// Represents a BAM file.
// Code is based heavily on bam.js, part of the Dalliance Genome Explorer,  (c) Thomas Down 2006-2001.

var igv = (function (igv) {

    var BAM_MAGIC = 21840194;
    var BAI_MAGIC = 21578050;
    var SECRET_DECODER = ['=', 'A', 'C', 'x', 'G', 'x', 'x', 'x', 'T', 'x', 'x', 'x', 'x', 'x', 'x', 'N'];
    var SECRET_DECODER_INT = new Uint8Array(SECRET_DECODER.map(function(s) { return s.charCodeAt(0); }));
    var CIGAR_DECODER = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '?', '?', '?', '?', '?', '?', '?'];
    var READ_STRAND_FLAG = 0x10;
    var MATE_STRAND_FLAG = 0x20;
    var FIRST_OF_PAIR_FLAG = 0x40;
    var SECOND_OF_PAIR_FLAG = 0x80;
    var NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
    var READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    var DUPLICATE_READ_FLAG = 0x400;
    var SUPPLEMENTARY_FLAG = 0x800;

    const MAX_GZIP_BLOCK_SIZE = 65536;   //  APPARENTLY.  Where is this documented???
    const DEFAULT_SAMPLING_WINDOW_SIZE = 100;
    const DEFAULT_SAMPLING_DEPTH = 50;
    const MAXIMUM_SAMPLING_DEPTH = 2500;

    /**
     * Class for reading a bam file
     *
     * @param config
     * @constructor
     */
    igv.BamReader = function (config) {
        this.cache = new igv.PromiseCache();

        this.config = config;

        this.filter = config.filter || new igv.BamFilter();

        this.bamPath = config.url;
        // Todo - deal with Picard convention.  WHY DOES THERE HAVE TO BE 2?
        this.baiPath = config.indexURL || this.bamPath + ".bai"; // If there is an indexURL provided, use it!
        this.headPath = config.headURL || this.bamPath;


        this.samplingWindowSize = config.samplingWindowSize === undefined ? DEFAULT_SAMPLING_WINDOW_SIZE : config.samplingWindowSize;
        this.samplingDepth = config.samplingDepth === undefined ? DEFAULT_SAMPLING_DEPTH : config.samplingDepth;
        if(this.samplingDepth > MAXIMUM_SAMPLING_DEPTH) {
            igv.log("Warning: attempt to set sampling depth > maximum value of 2500");
            this.samplingDepth = MAXIMUM_SAMPLING_DEPTH;
        }

        if (config.viewAsPairs) {
            this.pairsSupported = true;
        }
        else {
            this.pairsSupported = config.pairsSupported === undefined ? true : config.pairsSupported;
        }

    };

    igv.BamReader.prototype.readAlignments = function (chr, bpStart, bpEnd) {

        var self = this;


        return new Promise(function (fulfill, reject) {


            getChrIndex(self).then(function (chrToIndex) {

                var chrId = chrToIndex[chr],

                    alignmentContainer = new igv.AlignmentContainer(chr, bpStart, bpEnd, self.samplingWindowSize, self.samplingDepth, self.pairsSupported);


                if (!self.source) {
                    var options = {
                        headers: self.config.headers,
                        withCredentials: self.config.withCredentials
                    };

                    self.source =
                        new asyncArray.decompressBgzf(
                            new asyncArray.cache(new asyncArray.openUrlBuffer(self.bamPath, options)));
                }

                if (chrId === undefined) {
                    fulfill(alignmentContainer);
                } else {

                    getIndex(self).then(function (bamIndex) {

                        var chunks = bamIndex.blocksForRange(chrId, bpStart, bpEnd),
                            promises = [];

                        if (!chunks) {
                            fulfill(null);
                            reject("Error reading bam index");
                            return;
                        }
                        if (chunks.length === 0) {
                            fulfill(alignmentContainer);
                            return;
                        }

                        chunks.forEach(function (c) {

                            promises.push(new Promise(function (fulfill, reject) {

                                var fetchMin = c.minv.block,
                                    fetchMax = c.maxv.block;

                                self.source.read(fetchMin, fetchMax - fetchMin).then(function (blocks) {
                                    var alignments = [];
                                    var data = new Uint8Array(0);
                                    var minOffset = c.minv.offset;
                                    var i;
                                    var block;
                                    var outputBuffer = new Uint8Array(65536 * 2);
                                    var ba;

                                    for (i = 0; i < blocks.length; i++) {
                                        block = blocks[i];

                                        var blockAlignments = [];

                                        if (minOffset == 0 && data.byteLength == 0 && block.alignments) {
                                            blockAlignments = block.alignments;
                                        } else {
                                            var startTime  = performance.now();
                                            if (data.byteLength > 0) {
                                                ba = outputBuffer.subarray(0, data.byteLength + block.uncompressedSize);

                                                ba.set(data, 0);
                                                block.decompress(ba.subarray(data.byteLength));
                                            } else {
                                                ba = outputBuffer.subarray(0, block.uncompressedSize);
                                                block.decompress(ba);
                                            }
                                            var offset = decodeBamRecords(blockAlignments, ba, minOffset, self.filter);
                                            igv.addStat('decodeTime', performance.now() - startTime);
                                            igv.addStat('alignmentCount', blockAlignments.length);
                                            igv.addStat('data', ba.byteLength);
                                            if (minOffset == 0 && data.byteLength == 0 && offset == ba.byteLength) {
                                                block.alignments = blockAlignments;
                                            }
                                            data = ba.slice(offset);
                                        }

                                        alignments.push(blockAlignments.filter(function(alignment) {
                                            var refID = alignment.refID;
                                            var pos = alignment.start;

                                            if(refID < 0) {
                                                return 0;   // unmapped reads
                                            }
                                            else if (refID > chrId || pos > bpEnd) {
                                                return 0;    // off right edge
                                            }
                                            else if (refID < chrId || alignment.start + alignment.lengthOnRef < bpStart) {
                                                return 0;   // to left of start
                                            }
                                            return 1;
                                        }));

                                        if (blockAlignments.length > 0 && blockAlignments[blockAlignments.length - 1].start > bpEnd) {
                                            break;
                                        }
                                        minOffset = 0;
                                    }

                                    fulfill(alignments);

                                }).catch(function (obj) {
                                    reject(obj);
                                });

                            }));
                        });


                        Promise.all(promises).then(function (chunks) {
                            var startTime = performance.now();
                            chunks.forEach(function(alignments) {
                                alignments.forEach(function(blockAlignments) {
                                    blockAlignments.forEach(function(alignment) {
                                        alignmentContainer.push(alignment);
                                    });
                                });
                            });
                            alignmentContainer.finish();
                            igv.addStat('alignmentContainerTime', performance.now() - startTime);
                            fulfill(alignmentContainer);
                        }).catch(function (obj) {
                            reject(obj);
                        });
                    }).catch(reject);
                }
            }).catch(reject);
        });

        function decodeBamRecords(alignments, ba, offset, filter) {

            var blockSize,
                blockEnd,
                alignment,
                refID,
                pos,
                bmn,
                bin,
                mq,
                nl,
                flag_nc,
                flag,
                nc,
                lseq,
                mateRefID,
                matePos,
                j,
                k,
                p,
                seqBytes;

            while (offset + 4 <= ba.length) {

                blockSize = readInt(ba, offset);
                blockEnd = offset + blockSize + 4;

                if (blockEnd > ba.length) {
                    break;
                }

                alignment = new igv.BamAlignment();

                refID = readInt(ba, offset + 4);
                pos = readInt(ba, offset + 8);

                bmn = readInt(ba, offset + 12);
                bin = (bmn & 0xffff0000) >> 16;
                mq = (bmn & 0xff00) >> 8;
                nl = bmn & 0xff;

                flag_nc = readInt(ba, offset + 16);
                flag = (flag_nc & 0xffff0000) >> 16;
                nc = flag_nc & 0xffff;

                alignment.flags = flag;
                alignment.strand = !(flag & READ_STRAND_FLAG);

                lseq = readInt(ba, offset + 20);

                mateRefID = readInt(ba, offset + 24);
                matePos = readInt(ba, offset + 28);
                alignment.fragmentLength = readInt(ba, offset + 32);

                alignment.refID = refID;
                alignment.start = pos;
                alignment.mq = mq;
                alignment.chr = self.indexToChr[refID];

                if (!filter.pass(alignment)) {
                    offset = blockEnd;
                    continue;
                }

                alignment.readName = String.fromCharCode.apply(null, ba.subarray(offset + 36, offset + 36 + nl - 1));
                p = offset + 36 + nl;

                cigarOffset = p;
                p += nc * 4;

                seqBytes = (lseq + 1) >> 1;
                alignment.seq = new SequenceArray(ba.slice(p, p + seqBytes), lseq);
                p += seqBytes;

                if (lseq === 1 && String.fromCharCode(ba[p + j] + 33) === "*") {
                    // TODO == how to represent this?
                }
                else {
                    alignment.qual = ba.slice(p, p + lseq);
                }
                p += lseq;
                makeBlocks(alignment, ba, cigarOffset, nc);

                if (mateRefID >= 0) {
                    alignment.mate = {
                        chr: self.indexToChr[mateRefID],
                        position: matePos,
                        strand: !(flag & MATE_STRAND_FLAG)
                    };
                }


                alignment.tagBA = ba.slice(p, blockEnd);  // decode these on demand
                alignments.push(alignment);

                offset = blockEnd;
            }
            return offset;
            // Exits via top of loop.
        }

        /**
         * Split the alignment record into blocks as specified in the cigarArray.  Each aligned block contains
         * its portion of the read sequence and base quality strings.  A read sequence or base quality string
         * of "*" indicates the value is not recorded.  In all other cases the length of the block sequence (block.seq)
         * and quality string (block.qual) must == the block length.
         *
         * NOTE: Insertions are not yet treated // TODO
         *
         * @param record
         * @param cigarArray
         * @returns array of blocks
         */
        function makeBlocks(alignment, ba, p, len) {

            var blocks = [],
                insertions,
                seqOffset = 0,
                pos = alignment.start,
                seq = alignment.seq,
                qual = alignment.qual,
                gapType,
                cigar = '';

            for (var i = 0; i < len; i++) {
                var cigop = readInt(ba, p);
                var opLen = (cigop >> 4);
                var opLtr = CIGAR_DECODER[cigop & 0xf];

                cigar = cigar + opLen + opLtr;
                p += 4;

                switch (opLtr) {
                    case 'H' :
                        break; // ignore hard clips
                    case 'P' :
                        break; // ignore pads
                    case 'S' :
                        seqOffset += opLen;
                        gapType = 'S';
                        break; // soft clip read bases
                    case 'N' :
                        pos += opLen;
                        gapType = 'N';
                        break;  // reference skip
                    case 'D' :
                        pos += opLen;
                        gapType = 'D';
                        break;
                    case 'I' :
                        if (insertions === undefined) insertions = [];
                        insertions.push(makeBlock());
                        seqOffset += opLen;
                        break;
                    case 'M' :
                    case 'EQ' :
                    case '=' :
                    case 'X' :

                        blocks.push(makeBlock(gapType));
                        seqOffset += opLen;
                        pos += opLen;

                        break;

                    default :
                        console.log("Error processing cigar element: " + opLen + opLtr);
                }
            }
            alignment.cigar = cigar;
            alignment.blocks = blocks;
            alignment.insertions = insertions;
            alignment.lengthOnRef = pos - alignment.start;

            function makeBlock(gapType) {
                var blockSeq = seq;
                var blockQual = qual;

                if (seqOffset !== 0 || opLen !== seq.length) {
                    blockSeq = seq.subarray(seqOffset, seqOffset + opLen);
                    blockQual = qual && qual.subarray(seqOffset, seqOffset + opLen);
                }
                return new AlignmentBlock(pos, opLen, blockSeq, blockQual, gapType)
            }
        }
    }

    function AlignmentBlock(pos, len, seq, qual, gapType) {
        this.start = pos;
        this.len = len;
        this.seq = seq;
        this.qual = qual;
        this.gapType = gapType;
    }

    function SequenceArray(data, length, offset) {
        this.data = data;
        this.length = length;
        this.offset = offset || 0;
    }

    SequenceArray.prototype.charCodeAt = function(index) {
        var i = index + this.offset;
        var sb = this.data[i >> 1];

        if (i % 2 == 0) {
            return SECRET_DECODER_INT[(sb & 0xf0) >> 4]
        } else {
            return SECRET_DECODER_INT[(sb & 0x0f)]
        }
    };

    SequenceArray.prototype.charAt = function(index) {
        return String.fromCharCode(this.charCodeAt(index));
    };

    SequenceArray.prototype.toString = function() {
        var charCodes = new Array(this.length);
        var i;

        for (i = 0; i < this.length; i++) {
            charCodes[i] = this.charCodeAt(i);
        }
        return String.fromCharCode.apply(null, charCodes);
    };

    SequenceArray.prototype.subarray = function(start, end) {
        return new SequenceArray(this.data, end - start, this.offset + start);
    };

    igv.BamReader.prototype.readHeader = function () {

        var self = this;

        return new Promise(function (fulfill, reject) {

            getIndex(self).then(function (index) {

                var len = (index.firstAlignmentBlock == null ? 0 : index.firstAlignmentBlock) +
                    MAX_GZIP_BLOCK_SIZE; // Insure we get the complete compressed block containing the header

                igvxhr.loadArrayBuffer(self.bamPath,
                    {
                        headers: self.config.headers,

                        range: {start: 0, size: len},

                        withCredentials: self.config.withCredentials
                    }).then(function (compressedBuffer) {

                    var unc = igv.unbgzf(compressedBuffer, compressedBuffer.length),
                        uncba = new Uint8Array(unc),
                        magic = readInt(uncba, 0),
                        samHeaderLen = readInt(uncba, 4),
                        samHeader = '',
                        genome = igv.browser ? igv.browser.genome : null;

                    for (var i = 0; i < samHeaderLen; ++i) {
                        samHeader += String.fromCharCode(uncba[i + 8]);
                    }

                    var nRef = readInt(uncba, samHeaderLen + 8);
                    var p = samHeaderLen + 12;

                    self.chrToIndex = {};
                    self.indexToChr = [];
                    for (var i = 0; i < nRef; ++i) {
                        var lName = readInt(uncba, p);
                        var name = '';
                        for (var j = 0; j < lName - 1; ++j) {
                            name += String.fromCharCode(uncba[p + 4 + j]);
                        }
                        var lRef = readInt(uncba, p + lName + 4);
                        //dlog(name + ': ' + lRef);

                        if (genome && genome.getChromosomeName) {
                            name = genome.getChromosomeName(name);
                        }

                        self.chrToIndex[name] = i;
                        self.indexToChr.push(name);

                        p = p + 8 + lName;
                    }

                    fulfill();

                }).catch(reject);
            }).catch(reject);
        });
    }

//
    function getIndex(bam) {

        return new Promise(function (fulfill, reject) {

            if (bam.index) {
                fulfill(bam.index);
            }
            else {
                igv.loadBamIndex(bam.baiPath, bam.config).then(function (index) {
                    bam.index = index;

                    fulfill(bam.index);
                }).catch(reject);
            }
        });
    }


    function getChrIndex(bam) {

        return new Promise(function (fulfill, reject) {

            if (bam.chrToIndex) {
                fulfill(bam.chrToIndex);
            }
            else {
                bam.readHeader().then(function () {
                    fulfill(bam.chrToIndex);
                }).catch(reject);
            }
        });
    }

    function readInt(ba, offset) {
        return (ba[offset + 3] << 24) | (ba[offset + 2] << 16) | (ba[offset + 1] << 8) | (ba[offset]);
    }

    function readShort(ba, offset) {
        return (ba[offset + 1] << 8) | (ba[offset]);
    }

    return igv;

})
(igv || {});


