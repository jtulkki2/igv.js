/*
    Function for performing asynchronous reads
 */

var asyncArray = (function (asyncArray) {
    asyncArray.openUrlBuffer = function(url, options) {
        this.read = function(offset, length) {
            if (length == 0) {
                return Promise.resolve(new Uint8Array(0));
            }
            var blockOptions = Object.assign({}, options);

            blockOptions.range = {start: offset, size: length};
            return igvxhr.loadArrayBuffer(url, blockOptions);
        }
    }

    asyncArray.openUrlString = function(url, options) {
        this.read = function(offset, length) {
            if (length == 0) {
                return Promise.resolve("");
            }
            var blockOptions = Object.assign({}, options);

            blockOptions.range = {start: offset, size: length};
            return igvxhr.loadString(url, blockOptions);
        }
    }

    asyncArray.cache = function(array) {
        var cachedPages = [];

        function Page(start, end, data) {
            this.start = start;
            this.end = end;
            this.data = data;
        }

        function byPosition(a, b) {
            return a.start - b.start;
        }

        function getOverlappingPages(start, end) {
            return cachedPages.filter(function(page) {
                return page.start < end && page.end > start;
            });
        }

        function fetchPage(page) {
            page.data = array.read(page.start, page.end - page.start);
            cachedPages.push(page);
            return page.data;
        }

        function copyBytes(dest, src, offset) {
            var start = Math.max(offset, 0);
            var end = Math.min(offset + src.length, dest.length);

            if (start < end) {
                dest.set(src.subarray(start - offset, end - offset), start);
            }
        }

        function readDataFromPages(pages, start, length) {
            var data = new ArrayBuffer(length);
            var view = new Uint8Array(data);

            var promises = pages.map(function(page) {
                return page.data.then(function(data) {
                    copyBytes(view, new Uint8Array(data), page.start - start);
                });
            });

            return Promise.all(promises).then(function() {
                return data;
            });
        }

        this.read = function(offset, length) {
            var start = offset;
            var end = offset + length;
            var pages = getOverlappingPages(start, end).sort(byPosition);

            if (pages.length == 0) {
                return fetchPage(new Page(start, end));
            }
            var emptyStart = start;
            var emptyPages = [];
            pages.forEach(function(page) {
                if (page.start > emptyStart) {
                    emptyPages.push(new Page(emptyStart, page.start));
                }
                emptyStart = page.end;
            });
            if (emptyStart < end) {
                emptyPages.push(new Page(emptyStart, end));
            }

            emptyPages.forEach(fetchPage);
            pages = pages.concat(emptyPages);

            return readDataFromPages(pages, start, length);
        }
    }

    asyncArray.combineReads = function(array) {
        var chunkOffset;
        var chunkLength;
        var chunkPromise;
        var timeout;
        var timeoutFunc;

        function newChunk(offset, length) {
            if (timeout) {
                window.clearTimeout(timeout);
                timeoutFunc();
            }
            chunkOffset = offset;
            chunkLength = length;
            chunkPromise = new Promise(chunkExecutor);
        }

        function chunkExecutor(resolve, reject) {
            timeoutFunc = function() {
                var offset = chunkOffset;
                var length = chunkLength;

                timeout = 0;
                array.read(offset, length).then(function(data) {
                    chunkPromise = null;
                    resolve({
                        data: data,
                        offset: offset,
                        length: length
                    })
                }, function(error) {
                    chunkPromise = null;
                    reject(error)
                });
            }
            timeout = window.setTimeout(timeoutFunc);
        }

        this.read = function(offset, length) {
            if (chunkPromise && offset == chunkOffset + chunkLength) {
                chunkLength += length;
            } else {
                newChunk(offset, length);
            }

            return chunkPromise.then(function(chunk) {
                if (offset === chunk.offset && length === chunk.length) {
                    return chunk.data;
                }
                return chunk.data.slice(offset - chunk.offset, offset - chunk.offset + length);
            });
        }
    }

    asyncArray.decompressBgzf = function(array) {
        const MAX_GZIP_BLOCK_SIZE = 65536;
        var cachedBlocks = [];
        var cachedBlocksMap = {};

        function getBlocksInRange(start, end) {
            return cachedBlocks.filter(function(block) {
                return block.start >= start && block.start <= end;
            });
        }

        function addCachedBlock(block) {
            cachedBlocks.push(block);
            cachedBlocksMap[block.start] = block;
        }

        function byPosition(a, b) {
            return a.start - b.start;
        }

        function fetchBlocks(range) {
            var offset = range.start;
            var length = range.end - range.start;

            return array.read(offset, length + MAX_GZIP_BLOCK_SIZE).then(function(data) {
                var blocks = igv.unbgzfSplitBlocks(data);

                var startTime  = performance.now();
                blocks.forEach(function(block) {
                    if (block.start <= length && !cachedBlocksMap[block.start + offset]) {
                        igv.addStat('uncompressLength', block.end - block.start);
                        addCachedBlock({
                            start: block.start + offset,
                            end: block.end + offset,
                            data: jszlib_inflate_buffer_monolithic(data, block.start + 18, block.end - block.start - 18, [0], block.inputLength)
                        });
                    }
                });
                igv.addStat('uncompressTime', performance.now() - startTime);
                return blocks;
            });
        }

        function readDataFromBlocks(blocks) {
            var length = blocks.reduce(function(a, b) { return a + b.data.byteLength; }, 0);
            var data = new ArrayBuffer(length);
            var view = new Uint8Array(data);
            var offset = 0;

            blocks.forEach(function(block) {
                view.set(new Uint8Array(block.data), offset);
                offset += block.data.byteLength;
            });

            return data;
        }

        this.read = function(offset, length) {
            var start = offset;
            var end = offset + length;
            var blocks = getBlocksInRange(start, end).sort(byPosition);
            var emptyStart = start;
            var emptyRanges = [];

            blocks.forEach(function(block) {
                if (block.start > emptyStart) {
                    emptyRanges.push({start: emptyStart, end: block.start});
                }
                emptyStart = block.end;
            });
            if (emptyStart < end) {
                emptyRanges.push({start: emptyStart, end: end});
            }

            return Promise.all(emptyRanges.map(fetchBlocks)).then(function() {
                return getBlocksInRange(start, end).sort(byPosition);
//                return readDataFromBlocks(blocks, offset, length);
            });
        }
    }

    return asyncArray;
})
(asyncArray || {});
