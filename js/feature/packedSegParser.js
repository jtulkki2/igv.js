/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Broad Institute
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

/**
 *  Define parser for seg files  (.bed, .gff, .vcf, etc).  A parser should implement 2 methods
 *
 *     parseHeader(data) - return an object representing a header.  Details are format specific
 *
 *     parseFeatures(data) - return a list of features
 *
 */


var igv = (function (igv) {
    "use strict";

    var maxFeatureCount = Number.MAX_VALUE,    // For future use,  controls downsampling
        chrColumn = 0,
        startColumn = 1,
        endColumn = 2;


    igv.PackedSegParser = function () {
   }

    igv.PackedSegParser.prototype.parseHeader = function (data) {

        var lines = data.splitLines(),
            len = lines.length,
            line,
            i,
            tokens;

        for (i = 0; i < len; i++) {
            line = lines[i];
            if (line.startsWith("#")) {
                continue;
            }
            else {
                tokens = line.split("\t");
                this.header = {headings: tokens, lineCount: i + 1, samples: tokens.slice(3)};
                return this.header;
                break;
            }
        }

        return this.header;
    }


    igv.PackedSegParser.prototype.parseFeatures = function (data) {

        var lines = data ? data.splitLines() : [] ,
            len = lines.length,
            tokens, allFeatures = [], line, i, j, samples;

        if (!this.header) {
            this.header = this.parseHeader(data);
        }

        samples = this.header.samples;

        for (i = this.header.lineCount; i < len; i++) {

            line = lines[i];

            tokens = lines[i].split("\t");

            if (tokens.length == 3 + samples.length) {
                for (j = 0; j < samples.length; j++) {
                    allFeatures.push({
                        sample: samples[j],
                        chr: fixChromosome(tokens[chrColumn]),
                        start: parseInt(tokens[startColumn]),
                        end: parseInt(tokens[endColumn]),
                        value: parseValue(tokens[3 + j])
                    });
                }
            }
        }

        return allFeatures;

        // HOTFIX to fix broken SEG-files (added by Jouni 2018/2/13)
        function fixChromosome(chr) {
            if (chr == 23)
                return 'X';
            if (chr == 24)
                return 'Y';
            return chr;

        }

        function parseValue(str) {
            var value = parseFloat(str);

            return isNaN(value) ? null : value;
        }
    }


    return igv;
})
(igv || {});