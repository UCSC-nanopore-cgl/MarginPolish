//
// Created by Benedict Paten on 3/14/20.
//

#include "margin.h"

RleString * rleString_construct(char *str) {
    RleString *rleString = st_calloc(1, sizeof(RleString));

    rleString->nonRleLength = strlen(str);

    // Calc length of rle'd str
    for (uint64_t i = 0; i < rleString->nonRleLength; i++) {
        if (i + 1 == rleString->nonRleLength || str[i] != str[i + 1]) {
            rleString->length++;
        }
    }

    // Allocate
    rleString->rleString = st_calloc(rleString->length + 1, sizeof(char));
    rleString->repeatCounts = st_calloc(rleString->length, sizeof(uint64_t));

    // Fill out
    uint64_t j = 0, k = 1;
    for (uint64_t i = 0; i < rleString->nonRleLength; i++) {
        if (i + 1 == rleString->nonRleLength || str[i] != str[i + 1]) {
            rleString->rleString[j] = str[i];
            rleString->repeatCounts[j++] = k;
            k = 1;
        } else {
            k++;
        }
    }
    rleString->rleString[j] = '\0';
    assert(j == rleString->length);

    return rleString;
}


RleString *rleString_constructPreComputed(char *rleChars, const uint8_t *rleCounts) {
    RleString *rleString = st_calloc(1, sizeof(RleString));

    rleString->length = strlen(rleChars);
    rleString->nonRleLength = 0;
    for (int64_t i = 0; i < rleString->length; i++) {
        assert(rleCounts[i] >= 1);
        rleString->nonRleLength += rleCounts[i];
    }

    // Allocate
    rleString->rleString = stString_copy(rleChars);
    rleString->repeatCounts = st_calloc(rleString->length, sizeof(int64_t));

    // Fill out
    for (int64_t r = 0; r < rleString->length; r++) {
        // counts
        rleString->repeatCounts[r] = rleCounts[r];
    }

    return rleString;
}

RleString *rleString_construct_no_rle(char *string) {
    RleString *rleString = st_calloc(1, sizeof(RleString));

    rleString->nonRleLength = strlen(string);
    rleString->length = rleString->nonRleLength;

    // Allocate
    rleString->rleString = stString_copy(string);
    rleString->repeatCounts = st_calloc(rleString->length, sizeof(uint64_t));

    // Fill out repeat counts
    for (uint64_t i = 0; i < rleString->length; i++) {
        rleString->repeatCounts[i] = 1;
    }

    return rleString;
}

RleString *rleString_copySubstring(RleString *rleString, uint64_t start, uint64_t length) {
    RleString *rleSubstring = st_calloc(1, sizeof(RleString));

    assert(start + length <= rleString->length);

    // Set length of substring
    rleSubstring->length = length;

    // Copy character substring
    rleSubstring->rleString = stString_getSubString(rleString->rleString, start, length);

    // Copy repeat count substring and calculate non-rle length
    rleSubstring->nonRleLength = 0;
    rleSubstring->repeatCounts = st_calloc(length, sizeof(uint64_t));
    for (uint64_t i = 0; i < rleSubstring->length; i++) {
        rleSubstring->repeatCounts[i] = rleString->repeatCounts[i + start];
        rleSubstring->nonRleLength += rleSubstring->repeatCounts[i];
    }

    return rleSubstring;
}

void rleString_print(RleString *rleString, FILE *f) {
    fprintf(f, "%s -- ", rleString->rleString);
    for (int64_t i = 0; i < rleString->length; i++) {
        fprintf(f, "%" PRIi64 " ", rleString->repeatCounts[i]);
    }
}

RleString *rleString_copy(RleString *rleString) {
    return rleString_copySubstring(rleString, 0, rleString->length);
}

bool rleString_eq(RleString *r1, RleString *r2) {
    // If rle length or expanded lengths are not the same, then return false.
    if (r1->length != r2->length || r1->nonRleLength != r2->nonRleLength) {
        return 0;
    }
    // Check bases and repeat counts for equality
    for (int64_t i = 0; i < r1->length; i++) {
        if (r1->rleString[i] != r2->rleString[i] ||
            r1->repeatCounts[i] != r2->repeatCounts[i]) {
            return 0;
        }
    }
    return 1;
}

void rleString_destruct(RleString *rleString) {
    free(rleString->rleString);
    free(rleString->repeatCounts);
    free(rleString);
}

char *expandChar(char c, uint64_t repeatCount) {
    char *s = st_malloc(sizeof(char) * (repeatCount + 1));
    for (int64_t j = 0; j < repeatCount; j++) {
        s[j] = c;
    }
    s[repeatCount] = '\0';
    return s;
}

char *rleString_expand(RleString *rleString) {
    char *s = st_calloc(rleString->nonRleLength + 1, sizeof(char));
    int64_t j = 0;
    for (int64_t i = 0; i < rleString->length; i++) {
        for (int64_t k = 0; k < rleString->repeatCounts[i]; k++) {
            s[j++] = rleString->rleString[i];
        }
    }
    s[rleString->nonRleLength] = '\0';
    return s;
}

void rleString_rotateString(RleString *str, int64_t rotationLength, bool mergeEnds) {
    char rotatedString[str->length];
    uint64_t rotatedRepeatCounts[str->length];
    for (int64_t i = 0; i < str->length; i++) {
        rotatedString[(i + rotationLength) % str->length] = str->rleString[i];
        rotatedRepeatCounts[(i + rotationLength) % str->length] = str->repeatCounts[i];
    }
    int64_t j = 0;
    for (int64_t i = 0; i < str->length; i++) {
        if (!mergeEnds || i == 0 || rotatedString[i] != rotatedString[i - 1]) {
            str->rleString[j] = rotatedString[i];
            str->repeatCounts[j++] = rotatedRepeatCounts[i];
        } else {
            str->repeatCounts[j - 1] += rotatedRepeatCounts[i];
        }
        //str->rleString[i] = rotatedString[i];
        //str->repeatCounts[i] = rotatedRepeatCounts[i];
    }
    str->length = j;
}

uint8_t *rleString_rleQualities(RleString *rleString, const uint8_t *qualities) {
    // calculate read qualities (if set)
    //TODO unit test this
    uint8_t *rleQualities = st_calloc(rleString->length, sizeof(uint8_t));
    uint64_t rawPos = 0;
    for (uint64_t rlePos = 0; rlePos < rleString->length; rlePos++) {
        uint8_t min = UINT8_MAX;
        uint8_t max = 0;
        int64_t mean = 0;
        for (uint64_t repeatIdx = 0; repeatIdx < rleString->repeatCounts[rlePos]; repeatIdx++) {
            uint8_t q = qualities[rawPos++];
            min = (q < min ? q : min);
            max = (q > max ? q : max);
            mean += q;
        }
        mean = mean / rleString->repeatCounts[rlePos];
        assert(mean <= UINT8_MAX);
        // pick your favorite metric
        //r->qualities[rlePos] = min;
        //r->qualities[rlePos] = max;
        rleQualities[rlePos] = (uint8_t) mean;
    }
    assert(rawPos == rleString->nonRleLength);
    return rleQualities;
}

uint64_t *rleString_getNonRleToRleCoordinateMap(RleString *rleString) {
    uint64_t *nonRleToRleCoordinateMap = st_malloc(sizeof(uint64_t) * rleString->nonRleLength);

    uint64_t j = 0;
    for (uint64_t i = 0; i < rleString->length; i++) {
        for (uint64_t k = 0; k < rleString->repeatCounts[i]; k++) {
            nonRleToRleCoordinateMap[j++] = i;
        }
    }
    assert(j == rleString->nonRleLength);

    return nonRleToRleCoordinateMap;
}

uint64_t *rleString_getRleToNonRleCoordinateMap(RleString *rleString) {
    uint64_t *rleToNonRleCoordinateMap = st_malloc(sizeof(uint64_t) * rleString->length);

    uint64_t j = 0;
    for (uint64_t i = 0; i < rleString->length; i++) {
        rleToNonRleCoordinateMap[i] = j;
        j += rleString->repeatCounts[i];
    }
    assert(j == rleString->nonRleLength);

    return rleToNonRleCoordinateMap;
}

stList *runLengthEncodeAlignment(stList *alignment,
                                 const uint64_t *seqXNonRleToRleCoordinateMap,
                                 const uint64_t *seqYNonRleToRleCoordinateMap) {
    stList *rleAlignment = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    int64_t x = -1, y = -1;
    for (int64_t i = 0; i < stList_length(alignment); i++) {
        stIntTuple *alignedPair = stList_get(alignment, i);

        int64_t x2 = seqXNonRleToRleCoordinateMap[stIntTuple_get(alignedPair, 0)];
        int64_t y2 = seqYNonRleToRleCoordinateMap[stIntTuple_get(alignedPair, 1)];

        if (x2 > x && y2 > y) {
            stList_append(rleAlignment, stIntTuple_construct3(x2, y2, stIntTuple_get(alignedPair, 2)));
            x = x2;
            y = y2;
        }
    }

    return rleAlignment;
}