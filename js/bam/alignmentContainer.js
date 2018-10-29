/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
 * Author: Jim Robinson
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

var igv = (function (igv) {
    

    function canBePaired(alignment) {
        return alignment.isPaired() &&
            alignment.mate &&
            alignment.isMateMapped() &&
            alignment.chr === alignment.mate.chr &&
            (alignment.isFirstOfPair() || alignment.isSecondOfPair()) && !(alignment.isSecondary() || alignment.isSupplementary());
    }


    igv.AlignmentContainer = function (chr, start, end, samplingWindowSize, samplingDepth, pairsSupported) {

        this.chr = chr;
        this.start = start;
        this.end = end;
        this.length = (end - start);

        this.coverageMap = new CoverageMap(chr, start, end);
        this.alignments = [];
        this.downsampledIntervals = [];

        this.samplingWindowSize = samplingWindowSize === undefined ? 100 : samplingWindowSize;
        this.samplingDepth = samplingDepth === undefined ? 50 : samplingDepth;

        this.pairsSupported = pairsSupported;
        this.paired = false;  // false until proven otherwise
        this.pairsCache = {};  // working cache of paired alignments by read name

        this.downsampledReads = new Set();

        this.currentBucket = new DownsampleBucket(0, this.samplingWindowSize, this);

        this.filter = function filter(alignment) {         // TODO -- pass this in
            return alignment.isMapped() && !alignment.isFailsVendorQualityCheck();
        }

    }

    igv.AlignmentContainer.prototype.push = function (alignment) {

        if (this.filter(alignment) === false) return;

        this.coverageMap.incCounts(alignment);   // Count coverage before any downsampling

        if (this.pairsSupported && this.downsampledReads.has(alignment.readName)) {
            return;   // Mate already downsampled -- pairs are treated as a single alignment for downsampling
        }

        if (alignment.start >= this.currentBucket.end) {
            finishBucket.call(this);
            this.currentBucket = new DownsampleBucket(alignment.start, alignment.start + this.samplingWindowSize, this);
        }

        this.currentBucket.addAlignment(alignment);

    }

    igv.AlignmentContainer.prototype.forEach = function (callback) {
        this.alignments.forEach(callback);
    }

    igv.AlignmentContainer.prototype.finish = function () {

        if (this.currentBucket !== undefined) {
            finishBucket.call(this);
        }

        // Need to remove partial pairs whose mate was downsampled
        if(this.pairsSupported) {
            var tmp = [], ds = this.downsampledReads;

            this.alignments.forEach(function (a) {
                if (!ds.has(a.readName)) {
                    tmp.push(a);
                }
            })
            this.alignments = tmp;
        }

        this.alignments.sort(function (a, b) {
            return a.start - b.start
        });

        this.pairsCache = undefined;
        this.downsampledReads = undefined;
        this.coverageMap.finish();
    }

    igv.AlignmentContainer.prototype.contains = function (chr, start, end) {
        return this.chr == chr &&
            this.start <= start &&
            this.end >= end;
    }

    igv.AlignmentContainer.prototype.hasDownsampledIntervals = function () {
        return this.downsampledIntervals && this.downsampledIntervals.length > 0;
    }

    function finishBucket() {
        this.alignments = this.alignments.concat(this.currentBucket.alignments);
        if (this.currentBucket.downsampledCount > 0) {
            this.downsampledIntervals.push(new DownsampledInterval(
                this.currentBucket.start,
                this.currentBucket.end,
                this.currentBucket.downsampledCount));
        }
        this.paired = this.paired || this.currentBucket.paired;
    }

    function DownsampleBucket(start, end, alignmentContainer) {

        this.start = start;
        this.end = end;
        this.alignments = [];
        this.downsampledCount = 0;
        this.samplingDepth = alignmentContainer.samplingDepth;
        this.pairsSupported = alignmentContainer.pairsSupported;
        this.downsampledReads = alignmentContainer.downsampledReads;
        this.pairsCache = alignmentContainer.pairsCache;
    }

    DownsampleBucket.prototype.addAlignment = function (alignment) {

        var samplingProb, idx, replacedAlignment, pairedAlignment;

        if (this.pairsSupported && canBePaired(alignment) && (pairedAlignment = this.pairsCache[alignment.readName])) {
            // Not subject to downsampling, just update the existing alignment
            pairedAlignment.setSecondAlignment(alignment);
            this.pairsCache[alignment.readName] = undefined;   // Don't need to track this anymore. NOTE: Don't "delete", causes runtime performance issues
        } else if (this.alignments.length < this.samplingDepth) {

            if (this.pairsSupported && canBePaired(alignment)) {
                // First alignment in a pair
                pairedAlignment = new igv.PairedAlignment(alignment);
                this.paired = true;
                this.pairsCache[alignment.readName] = pairedAlignment;
                this.alignments.push(pairedAlignment);
            }
            else {
                this.alignments.push(alignment);
            }

        } else {

            samplingProb = this.samplingDepth / (this.samplingDepth + this.downsampledCount + 1);

            if (Math.random() < samplingProb) {

                idx = Math.floor(Math.random() * (this.alignments.length - 1));
                replacedAlignment = this.alignments[idx];   // To be replaced

                if (this.pairsSupported && canBePaired(alignment)) {

                    if(this.pairsCache[replacedAlignment.readName] !== undefined) {
                        this.pairsCache[replacedAlignment.readName] = undefined;
                    }

                    pairedAlignment = new igv.PairedAlignment(alignment);
                    this.paired = true;
                    this.pairsCache[alignment.readName] = pairedAlignment;
                    this.alignments[idx] = pairedAlignment;

                }
                else {
                    this.alignments[idx] = alignment;
                }
                this.downsampledReads.add(replacedAlignment.readName);

            }
            else {
                this.downsampledReads.add(alignment.readName);
            }

            this.downsampledCount++;
        }

    }

    var baseIdx = {
        A: 1,
        C: 2,
        T: 3,
        G: 4,
        N: 5
    };
    var baseCharCodeIdx = toCharCodeIdx(baseIdx);

    function toCharCodeIdx(stringIdx) {
        var idx = new Uint8Array(256);

        Object.keys(stringIdx).forEach(function(key) {
            idx[key.charCodeAt(0)] = stringIdx[key];
        });

        return idx;
    }

    // TODO -- refactor this to use an object, rather than an array,  if end-start is > some threshold
    function CoverageMap(chr, start, end) {

        var i;

        this.chr = chr;
        this.bpStart = start;
        this.length = (end - start);

        this.coverage = new Array(this.length);
        this.posByBase = new Array(6);
        this.negByBase = new Array(6);
        this.qualByBase = new Array(6);

        for (i = 0; i < 6; i++) {
            this.posByBase[i] = new Int32Array(this.length);
            this.negByBase[i] = new Int32Array(this.length);
            this.qualByBase[i] = new Int32Array(this.length);
        }

        this.maximum = 0;

        this.threshold = 0.2;
        this.qualityWeight = true;
    }

    CoverageMap.prototype.finish = function() {
        var self = this;
        var i, coverage, total, qual;

        for (i = 0; i < self.coverage.length; i++) {
            coverage = self.coverage[i] = new Coverage();
            total = 0;
            qual = 0;
            Object.keys(baseIdx).forEach(function(base) {
                var idx = baseIdx[base];

                total += coverage["pos" + base] = self.posByBase[idx][i];
                total += coverage["neg" + base] = self.negByBase[idx][i];
                qual += coverage["qual" + base] = self.qualByBase[idx][i];
            });
            total += self.posByBase[0][i];
            total += self.negByBase[0][i];
            qual += self.qualByBase[0][i];
            coverage.total = total;
            coverage.qual = qual;
            self.maximum = Math.max(coverage.total, self.maximum);
        }
        // free up memory
        this.posByBase = undefined;
        this.negByBase = undefined;
        this.qualByBase = undefined;
    }

    CoverageMap.prototype.incCounts = function (alignment) {

        var self = this;
        var coverageByBase = alignment.strand ? self.posByBase : self.negByBase;
        var qualByBase = self.qualByBase;


        if (alignment.blocks === undefined) {

            incBlockCount(alignment);
        }
        else {
            alignment.blocks.forEach(function (block) {
                incBlockCount(block);
            });
        }

        function incBlockCount(block) {

            var base,
                i,
                j,
                seq = block.seq,
                qual = block.qual;

            for (i = block.start - self.bpStart, j = 0; j < block.len; i++, j++) {

                base = baseCharCodeIdx[seq.charCodeAt(j)];

                coverageByBase[base][i]++;
                qualByBase[base][i] += qual[j];

            }
        }
    }

    function Coverage() {
        this.posA = 0;
        this.negA = 0;

        this.posT = 0;
        this.negT = 0;

        this.posC = 0;
        this.negC = 0;
        this.posG = 0;

        this.negG = 0;

        this.posN = 0;
        this.negN = 0;

        this.pos = 0;
        this.neg = 0;

        this.qualA = 0;
        this.qualT = 0;
        this.qualC = 0;
        this.qualG = 0;
        this.qualN = 0;

        this.qual = 0;

        this.total = 0;
    }

    Coverage.prototype.isMismatch = function (refBase) {

        var myself = this,
            mismatchQualitySum,
            threshold = igv.CoverageMap.threshold * ((igv.CoverageMap.qualityWeight && this.qual) ? this.qual : this.total);

        mismatchQualitySum = 0;
        ["A", "T", "C", "G"].forEach(function (base) {

            if (base !== refBase) {
                mismatchQualitySum += ((igv.CoverageMap.qualityWeight && myself.qual) ? myself["qual" + base] : (myself["pos" + base] + myself["neg" + base]));
            }
        });

        return mismatchQualitySum >= threshold;

    };

    DownsampledInterval = function (start, end, counts) {
        this.start = start;
        this.end = end;
        this.counts = counts;
    }

    DownsampledInterval.prototype.popupData = function (genomicLocation) {
        return [
            {name: "start", value: this.start + 1},
            {name: "end", value: this.end},
            {name: "# downsampled:", value: this.counts}]
    }


    return igv;

})(igv || {});