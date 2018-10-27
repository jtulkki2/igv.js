// Represents a BAM index.
// Code is based heavily on bam.js, part of the Dalliance Genome Explorer,  (c) Thomas Down 2006-2001.

var igv = (function (igv) {


    const BAI_MAGIC = 21578050;
    const TABIX_MAGIC = 21578324;
    const MAX_HEADER_SIZE = 100000000;   // IF the header is larger than this we can't read it !
    const MAX_GZIP_BLOCK_SIZE = (1 << 16);


    /**
     * @param indexURL
     * @param config
     * @param tabix
     *
     * @returns a Promised for the bam or tabix index.  The fulfill function takes the index as an argument.
     */
    igv.loadBamIndex = function (indexURL, config, tabix) {

        return new Promise(function (fulfill, reject) {

            var genome = igv.browser ? igv.browser.genome : null;

            igvxhr.loadArrayBuffer(indexURL,
                {
                    headers: config.headers,
                    withCredentials: config.withCredentials
                }).then(function (arrayBuffer) {

                var indices = [],
                    magic, nbin, nintv, nref, parser,
                    blockMin = Number.MAX_VALUE,
                    blockMax = null,
                    binIndex, linearIndex, binNumber, cs, ce, b, i, ref, sequenceIndexMap;

                if (!arrayBuffer) {
                    fulfill(null);
                    return;
                }

                if (tabix) {
                    var inflate = new Zlib.Gunzip(new Uint8Array(arrayBuffer));
                    arrayBuffer = inflate.decompress().buffer;
                }

                parser = new igv.BinaryParser(new DataView(arrayBuffer));

                magic = parser.getInt();

                if (magic === BAI_MAGIC || (tabix && magic === TABIX_MAGIC)) {

                    nref = parser.getInt();


                    if (tabix) {
                        // Tabix header parameters aren't used, but they must be read to advance the pointer
                        var format = parser.getInt();
                        var col_seq = parser.getInt();
                        var col_beg = parser.getInt();
                        var col_end = parser.getInt();
                        var meta = parser.getInt();
                        var skip = parser.getInt();
                        var l_nm = parser.getInt();

                        sequenceIndexMap = {};
                        for (i = 0; i < nref; i++) {
                            var seq_name = parser.getString();

                            // Translate to "official" chr name.
                            if (genome) seq_name = genome.getChromosomeName(seq_name);

                            sequenceIndexMap[seq_name] = i;
                        }
                    }

                    var blockSizes = {};
                    var prevBlock = 0;
                    const MAX_GZIP_BLOCK_SIZE = (1 << 16);

                    for (ref = 0; ref < nref; ++ref) {

                        binIndex = {};
                        linearIndex = [];

                        nbin = parser.getInt();

                        for (b = 0; b < nbin; ++b) {

                            binNumber = parser.getInt();

                            if (binNumber == 37450) {
                                // This is a psuedo bin, not used but we have to consume the bytes
                                nchnk = parser.getInt(); // # of chunks for this bin
                                cs = parser.getVPointer();   // unmapped beg
                                ce = parser.getVPointer();   // unmapped end
                                var n_maped = parser.getLong();
                                var nUnmapped = parser.getLong();

                            }
                            else {
                                
                                binIndex[binNumber] = [];
                                var nchnk = parser.getInt(); // # of chunks for this bin

                                for (i = 0; i < nchnk; i++) {
                                    cs = parser.getVPointer();
                                    ce = parser.getVPointer();
                                    if (cs && ce) {
                                        if (cs.block < blockMin) {
                                            blockMin = cs.block;    // Block containing first alignment
                                        }
                                        if (!blockMax || ce.block > blockMax.block ||
                                            ce.block == blockMax.block && ce.offset > blockMax.offset) {
                                            blockMax = ce;
                                        }
                                        blockSizes[cs.block] = null;
                                        blockSizes[ce.block] = null;
                                        binIndex[binNumber].push([cs, ce]);
                                    }
                                }
                            }
                        }


                        nintv = parser.getInt();
                        for (i = 0; i < nintv; i++) {
                            cs = parser.getVPointer();
                            linearIndex.push([cs, null]);   // Might be null
                            if (cs) {
                                blockSizes[cs.block] = null;
                            }
                        }

                        Object.keys(binIndex).forEach(function(bin) {
                            var index = binIndex[bin];
                            var reg = bin2reg(bin);
                            var i, j, lin, ce;

                            for (i = 0; i < index.length; i++) {
                                ce = index[i][1];
                                lin = linearIndex[(reg.beg >> 14) - 1];
                                if (lin && (!lin[1] || lin[1].block > ce.block ||
                                    (lin[1].block == ce.block && lin[1].offset > ce.offset))) {
                                    lin[1] = ce;
                                }
                            }
                        });

                        if (linearIndex.length > 0) {
                            linearIndex[linearIndex.length - 1][1] = blockMax;
                            for (i = linearIndex.length - 2; i >= 0; i--) {
                                if (!linearIndex[i][1]) {
                                    linearIndex[i][1] = linearIndex[i + 1][1];
                                }
                            }
                        }

                        if (nbin > 0) {
                            indices[ref] = {
                                binIndex: binIndex,
                                linearIndex: linearIndex
                            }
                        }
                    }
                    var blocks = Object.keys(blockSizes)
                        .map(function(block) { return +block; })
                        .sort(function(a, b) { return a - b; });

                    for (i = 0; i < blocks.length; i++) {
                        if (blocks[i] > prevBlock + MAX_GZIP_BLOCK_SIZE) {
                            blockSizes[prevBlock] = MAX_GZIP_BLOCK_SIZE;
                        } else {
                            blockSizes[prevBlock] = blocks[i] - prevBlock;
                        }
                        prevBlock = blocks[i];
                    }


                    if (blockSizes) {
                        blockSizes[prevBlock] = MAX_GZIP_BLOCK_SIZE;
                    }

                } else {
                    throw new Error(indexURL + " is not a " + (tabix ? "tabix" : "bai") + " file");
                }
                fulfill(new igv.BamIndex(indices, blockMin, sequenceIndexMap, tabix, blockSizes, blocks));
            }).catch(reject);
        })
    }


    igv.BamIndex = function (indices, blockMin, sequenceIndexMap, tabix, blockSizes, blocks) {
        this.firstAlignmentBlock = blockMin;
        this.indices = indices;
        this.sequenceIndexMap = sequenceIndexMap;
        this.tabix = tabix;
        this.blockSizes = blockSizes;
        this.blocks = blocks;

    }

    /**
     * Fetch blocks for a particular genomic range.  This method is public so it can be unit-tested.
     *
     * @param refId  the sequence dictionary index of the chromosome
     * @param min  genomic start position
     * @param max  genomic end position
     * @param return an array of {minv: {filePointer, offset}, {maxv: {filePointer, offset}}
     */
    igv.BamIndex.prototype.blocksForRange = function (refId, min, max) {

        var bam = this,
            ba = bam.indices[refId],
            minLin,
            maxLin,
            i,
            intChunks = [];

        if (!ba) {
            return intChunks;
        }
        else {

            l = ba.linearIndex;
            minLin = Math.min(min >> 14, l.length - 1);
            maxLin = Math.min(max >> 14, l.length - 1);
            var chunk;

            intChunks.push({
                minv: l[minLin][0],
                maxv: l[maxLin][1]
            });
/*            for (i = minLin; i <= maxLin; ++i) {
                chunk = l[i];
                var cs = chunk[0];
                var ce = chunk[1];
                intChunks.push({
                    minv: cs,
                    maxv: ce
                });
            }
*/

            return intChunks;
        }

    };


    /**
     * Calculate the list of bins that may overlap with region [beg, end]
     *
     */
    function reg2bins(beg, end) {
        var i = 0, k, list = [];
        if (end >= 1 << 29)   end = 1 << 29;
        --end;
        list.push(0);
        for (k = 1 + (beg >> 26); k <= 1 + (end >> 26); ++k) list.push(k);
        for (k = 9 + (beg >> 23); k <= 9 + (end >> 23); ++k) list.push(k);
        for (k = 73 + (beg >> 20); k <= 73 + (end >> 20); ++k) list.push(k);
        for (k = 585 + (beg >> 17); k <= 585 + (end >> 17); ++k) list.push(k);
        for (k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); ++k) list.push(k);
        return list;
    }

    function bin2reg(bin) {
        var length = 1 << 29;
        var nbin = 1;

        while (length >= 1 << 14) {
            if (bin < nbin) {
                return {beg: bin * length, end: (bin + 1) * length};
            }
            bin -= nbin;
            nbin <<= 3;
            length >>= 3;
        }
        return undefined;
    }

    return igv;

})(igv || {});