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
                    blockMin = null,
                    blockMax = 0,
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

                    for (ref = 0; ref < nref; ++ref) {

                        binIndex = {};
                        linearIndex = [];

                        nbin = parser.getInt();

                        var binPosition = parser.position;

                        for (b = 0; b < nbin; ++b) {

                            binNumber = parser.getInt();

                            if (binNumber == 37450) {
                                parser.position += 4 + 8 * 4;
                            }
                            else {
                                binIndex[binNumber] = [];
                                var nchnk = parser.getInt(); // # of chunks for this bin

                                parser.position += nchnk * 8 * 2;
                            }
                        }

                        nintv = parser.getInt();
                        linearStartOffset = new Uint32Array(nintv);
                        linearStartBlock = new Float64Array(nintv);
                        linearEndBlock = new Float64Array(nintv);
                        for (i = 0; i < nintv; i++) {
                            linearStartOffset[i] = parser.getUShort();
                            linearStartBlock[i] = parser.getBlock();
                        }

                        var endPosition = parser.position;
                        parser.position = binPosition;

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
                                var reg = bin2reg(binNumber);
                                var j = (reg.beg >> 14) - 1;

                                for (i = 0; i < nchnk; i++) {
                                    var startOffset = parser.getUShort();
                                    var startBlock = parser.getBlock();
                                    var endOffset = parser.getUShort();
                                    var endBlock = parser.getBlock();

                                    if (blockMin == null || startBlock < blockMin) {
                                        blockMin = startBlock;    // Block containing first alignment
                                    }
                                    if (endBlock > blockMax) {
                                        blockMax = endBlock;
                                    }
                                    if (j >= 0 && j < linearEndBlock.length && (!linearEndBlock[j] || startBlock < linearEndBlock[j])) {
                                        linearEndBlock[j] = startBlock + 65536;
                                    }
                                }
                            }
                        }
                        parser.position = endPosition;

                        if (linearEndBlock.length > 0) {
                            linearEndBlock[linearEndBlock.length - 1] = blockMax;
                            for (i = linearEndBlock.length - 2; i >= 0; i--) {
                                if (!linearEndBlock[i]) {
                                    linearEndBlock[i] = linearEndBlock[i + 1];
                                }
                            }
                        }

                        if (nbin > 0) {
                            indices[ref] = {
                                linearStartOffset: linearStartOffset,
                                linearStartBlock: linearStartBlock,
                                linearEndBlock: linearEndBlock
                            }
                        }
                    }

                } else {
                    throw new Error(indexURL + " is not a " + (tabix ? "tabix" : "bai") + " file");
                }
                fulfill(new igv.BamIndex(indices, blockMin, sequenceIndexMap, tabix));
            }).catch(reject);
        })
    }


    igv.BamIndex = function (indices, blockMin, sequenceIndexMap, tabix) {
        this.firstAlignmentBlock = blockMin;
        this.indices = indices;
        this.sequenceIndexMap = sequenceIndexMap;
        this.tabix = tabix;
//        this.blockSizes = blockSizes;
//        this.blocks = blocks;

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

            var length = ba.linearStartBlock.length;
            minLin = Math.min(min >> 14, length - 1);
            maxLin = Math.min(max >> 14, length - 1);

            intChunks.push({
                minv: {
                    block: ba.linearStartBlock[minLin],
                    offset: ba.linearStartOffset[minLin]
                },
                maxv: {
                    block: ba.linearEndBlock[maxLin],
                    offset: 0
                }
            });

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