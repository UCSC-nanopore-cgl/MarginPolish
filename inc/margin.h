/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef ST_RP_HMM_H_
#define ST_RP_HMM_H_

#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <ctype.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "sonLib.h"
#include "hashTableC.h"
#include "pairwiseAligner.h"
#include "randomSequences.h"
#include "stateMachine.h"

#include <htslib/hts.h>
#include <htslib/sam.h>

# ifdef _OPENMP
#include <omp.h>
# endif


#define uint128_t __uint128_t

/*
 * More function documentation is in the .c files
 */

/*
 * Phasing structs
 */
typedef struct _stSite stSite;
typedef struct _stReference stReference;
typedef struct _stProfileSeq stProfileSeq;
typedef struct _stRPHmm stRPHmm;
typedef struct _stRPHmmParameters stRPHmmParameters;
typedef struct _stRPColumn stRPColumn;
typedef struct _stRPCell stRPCell;
typedef struct _stRPMergeColumn stRPMergeColumn;
typedef struct _stRPMergeCell stRPMergeCell;
typedef struct _stGenomeFragment stGenomeFragment;

/*
 * Polisher structs
 */
typedef struct _stReadHaplotypeSequence stReadHaplotypeSequence;
typedef struct hashtable stReadHaplotypePartitionTable;
typedef struct _repeatSubMatrix RepeatSubMatrix;
typedef struct _polishParams PolishParams;
typedef struct _Poa Poa;
typedef struct _poaNode PoaNode;
typedef struct _poaInsert PoaInsert;
typedef struct _poaDelete PoaDelete;
typedef struct _poaBaseObservation PoaBaseObservation;
typedef struct _rleString RleString;
typedef struct _refMsaView MsaView;
typedef struct _outputChunkers OutputChunkers;

/*
 * VCF structs
 */
typedef struct _vcfEntry VcfEntry;

/*
 * Combined params object
 */
typedef struct _params Params;

/*
 * Combined parameter object for phase, polish, view, etc.
 */

struct _params {
	PolishParams *polishParams;
	stRPHmmParameters *phaseParams;
};

Params *params_readParams(char *paramsFile);

void params_destruct(Params *params);

void params_printParameters(Params *params, FILE *fh);

/*
 * Overall coordination functions
 */

stList *filterReadsByCoverageDepth(stList *profileSeqs, stRPHmmParameters *params, stList *filteredProfileSeqs,
								   stList *discardedProfileSeqs);

stList *getRPHmms(stList *profileSeqs, stRPHmmParameters *params);

stList *getTilingPaths(stSortedSet *hmms);

stSet *getOverlappingComponents(stList *tilingPath1, stList *tilingPath2);

stList *mergeTwoTilingPaths(stList *tilingPath1, stList *tilingPath2);

stRPHmm *fuseTilingPath(stList *tilingPath);

/*
 * Math
 */
#define ST_MATH_LOG_ZERO -INFINITY
#define ST_MATH_LOG_ONE 0.0

double logAddP(double a, double b, bool maxNotSum);

/*
 * Strandedness
 */
#define POS_STRAND_IDX 1
#define NEG_STRAND_IDX 0

/*
 * Repeat Count
 */

#define MAXIMUM_REPEAT_LENGTH 51

#define ALLELE_LOG_PROB_BITS 8

/*
 * Binary partition stuff
 */

// The maximum read depth the model can support
#define MAX_READ_PARTITIONING_DEPTH 64

char *intToBinaryString(uint64_t i);

uint64_t makeAcceptMask(uint64_t depth);

uint64_t mergePartitionsOrMasks(uint64_t partition1, uint64_t partition2,
								uint64_t depthOfPartition1, uint64_t depthOfPartition2);

uint64_t maskPartition(uint64_t partition, uint64_t mask);

bool seqInHap1(uint64_t partition, int64_t seqIndex);

uint64_t invertPartition(uint64_t partition, uint64_t depth);

uint64_t flipAReadsPartition(uint64_t partition, uint64_t readIndex);


/*
 * Reference / site definition
 */

struct _stSite {
	uint64_t alleleNumber; // Number of alleles at the site in the reference
	uint64_t alleleOffset; // The index of the first allele in this site
	// in a sequence of all alleles in the reference, ordered first by site then
	// by order in the site.
	uint16_t *substitutionLogProbs; // Log probabilities of substitutions between the alleles
	uint16_t *allelePriorLogProbs; // Prior log-prob on alleles, allows upweighting of reference allele
};

uint16_t *stSite_getSubstitutionProb(stSite *s, int64_t from, int64_t to);

struct _stReference {
	char *referenceName;
	uint64_t length; // Number of sites
	uint64_t totalAlleles; // Total number of alleles across all sites
    stSite *sites;
};

void stReference_destruct(stReference *ref);

/*
 * _stProfileSeq
 * Struct for profile sequence
 */

#define PROFILE_PROB_SCALAR 30.0

struct _stProfileSeq {
    stReference *ref;
    char *readId;
    uint64_t refStart; // The first site in the reference
    uint64_t length; // Number of reference sites
    uint64_t alleleOffset; // The index of the first allele in this sequence
    // in a sequence of all alleles in the reference, ordered first by site then
    // by order in the site.

	// The log-probability of alleles, as specified by uint8_t
	// Each is expressed as an 8 bit unsigned int, with the value from 0 to -255
	uint8_t *profileProbs;
};

stProfileSeq *stProfileSeq_constructEmptyProfile(stReference *ref, char *readId,
												 int64_t referenceStart, int64_t length);

void stProfileSeq_destruct(stProfileSeq *seq);

void stProfileSeq_print(stProfileSeq *seq, FILE *fileHandle);

uint8_t *stProfileSeq_getProb(stProfileSeq *seq, uint64_t site, uint64_t allele);

/*
 * Emission probabilities methods
 */
double emissionLogProbability(stRPColumn *column, stRPCell *cell, uint64_t *bitCountVectors,
							  stReference *reference,
							  stRPHmmParameters *params);

void fillInPredictedGenome(stGenomeFragment *gF, uint64_t partition,
						   stRPColumn *column, stRPHmmParameters *params);

/*
 * Constituent functions tested and used to do bit twiddling
*/
int popcount64(uint64_t x);

uint64_t getLogProbOfAllele(uint64_t *bitCountVectors, uint64_t depth, uint64_t partition,
							uint64_t siteOffset, uint64_t allele);

uint64_t *calculateCountBitVectors(uint8_t **seqs, stReference *ref,
								   uint64_t firstSite, uint64_t length, uint64_t depth);

/*
 * _stRPHmmParameters
 * Struct for hmm parameters
 */
struct _stRPHmmParameters {
	/*
	 * Parameters used for the HMM computation
	 */
	bool maxNotSumTransitions;

	// Filters on the number of states in a column
	// Used to prune the hmm
	int64_t minPartitionsInAColumn;
	int64_t maxPartitionsInAColumn;
	double minPosteriorProbabilityForPartition;

	// MaxCoverageDepth is the maximum depth of profileSeqs to allow at any base.
	// If the coverage depth is higher than this then some profile seqs are randomly discarded.
	int64_t maxCoverageDepth;
	int64_t minReadCoverageToSupportPhasingBetweenHeterozygousSites;

	// Ensure symmetry in the HMM such that the inverted partition of each partition is included in the HMM
	bool includeInvertedPartitions;

	// Number of rounds of iterative refinement to attempt to improve the partition.
	int64_t roundsOfIterativeRefinement;

	// Flag used to determine if the ancestor substitution probabilities are used
	bool includeAncestorSubProb;

	// Whether to include reads in a partition
	int64_t minPhredScoreForHaplotypePartition;

	// Should homozygous variants be used during bubble finding with input VCF? (improves sequence quality, can confound phasing)
	bool includeHomozygousVCFEntries;

    // should we include all vcf entries, or just ["PASS", "pass", "."]
    bool onlyUsePassVCFEntries;

    // should we include all vcf entries, or just SNPs
    bool onlyUseSNPVCFEntries;

	// should we stitch with primary reads or all reads (including filtering via downsampling or other means)
	bool stitchWithPrimaryReadsOnly;

	// should we use sampling, where we keep..
	bool useVariantSelectionAdaptiveSampling;

	// ..all variants above this threshold
	double variantSelectionAdaptiveSamplingPrimaryThreshold;

	// ..and if we don't have (on average) a variant every this many bp
	int64_t variantSelectionAdaptiveSamplingDesiredBasepairsPerVariant;

	// ..then take variants sorted by qual until we hit that threshold. but not any variant below this threshold:
    double minVariantQuality;

    // for output vcf
    bool updateAllOutputVCFFormatFields;

    // likelihood at which we start a new phase set based on read concordance between haplotypes
    double phasesetMinBinomialReadSplitLikelihood;

    // ratio of discordant reads to all reads for phase set determination
    double phasesetMaxDiscordantRatio;

	// Number of iterations to search for bubbles (and remove bubbles with strand or read split below some threshold)
	int64_t bubbleFindingIterations;

	// Parameters for removing bubbles
	double bubbleMinBinomialStrandLikelihood;
	double bubbleMinBinomialReadSplitLikelihood;

};

stRPHmmParameters *stRPHmmParameters_construct();
stRPHmmParameters *stRPHmmParameters_copy(stRPHmmParameters *toCopy);
void stRPHmmParameters_destruct(stRPHmmParameters *params);

void stRPHmmParameters_printParameters(stRPHmmParameters *params, FILE *fH);

/*
 * _stRPHmm
 * Struct for read partitioning hmm
 */
struct _stRPHmm {
	stReference *ref;
	int64_t refStart; // First site in reference
	int64_t refLength; // Number of sites in the reference
	stList *profileSeqs; // List of stProfileSeq
	int64_t columnNumber; // Number of columns, excluding merge columns
	int64_t maxDepth;
	stRPColumn *firstColumn;
	stRPColumn *lastColumn;
	const stRPHmmParameters *parameters;
	//Forward/backward probability calculation things
	double forwardLogProb;
	double backwardLogProb;
};

stRPHmm *stRPHmm_construct(stProfileSeq *profileSeq, stRPHmmParameters *params);

void stRPHmm_destruct(stRPHmm *hmm, bool destructColumns);

void stRPHmm_destruct2(stRPHmm *hmm);

bool stRPHmm_overlapOnReference(stRPHmm *hmm1, stRPHmm *hmm2);

stRPHmm *stRPHmm_createCrossProductOfTwoAlignedHmm(stRPHmm *hmm1, stRPHmm *hmm2);

void stRPHmm_alignColumns(stRPHmm *hmm1, stRPHmm *hmm2);

stRPHmm *stRPHmm_fuse(stRPHmm *leftHmm, stRPHmm *rightHmm);

void stRPHmm_forwardBackward(stRPHmm *hmm);

void stRPHmm_prune(stRPHmm *hmm);

void stRPHmm_print(stRPHmm *hmm, FILE *fileHandle, bool includeColumns, bool includeCells);

stList *stRPHmm_forwardTraceBack(stRPHmm *hmm);

stSet *stRPHmm_partitionSequencesByStatePath(stRPHmm *hmm, stList *path, bool partition1);

int stRPHmm_cmpFn(const void *a, const void *b);

stRPHmm *stRPHmm_split(stRPHmm *hmm, int64_t splitPoint);

void stRPHmm_resetColumnNumberAndDepth(stRPHmm *hmm);

stList *stRPHMM_splitWherePhasingIsUncertain(stRPHmm *hmm);

void logHmm(stRPHmm *hmm, stGenomeFragment *gF);

/*
 * _stRPColumn
 * Column of read partitioning hmm
 */
struct _stRPColumn {
	int64_t refStart; // First site in the reference
	int64_t length; // Number of sites
	int64_t depth;
	stProfileSeq **seqHeaders;
	uint8_t **seqs;
	stRPCell *head;
	stRPMergeColumn *nColumn, *pColumn;
	double totalLogProb;
};

stRPColumn *stRPColumn_construct(int64_t refStart, int64_t length, int64_t depth,
								 stProfileSeq **seqHeaders, uint8_t **seqs);

void stRPColumn_destruct(stRPColumn *column);

void stRPColumn_print(stRPColumn *column, FILE *fileHandle, bool includeCells);

void stRPColumn_split(stRPColumn *column, int64_t firstHalfLength, stRPHmm *hmm);

stSet *stRPColumn_getSequencesInCommon(stRPColumn *column1, stRPColumn *column2);

stSet *stRPColumn_getColumnSequencesAsSet(stRPColumn *column);

/*
 * _stRPCell
 * State of read partitioning hmm
 */
struct _stRPCell {
	uint64_t partition;
	double forwardLogProb, backwardLogProb;
	stRPCell *nCell;
};

stRPCell *stRPCell_construct(int64_t partition);

void stRPCell_destruct(stRPCell *cell);

void stRPCell_print(stRPCell *cell, FILE *fileHandle);

double stRPCell_posteriorProb(stRPCell *cell, stRPColumn *column);

/*
 * _stRPMergeColumn
 * Merge column of read partitioning hmm
 */
struct _stRPMergeColumn {
	uint64_t maskFrom;
	uint64_t maskTo;
	stHash *mergeCellsFrom;
	stHash *mergeCellsTo;
	stRPColumn *nColumn, *pColumn;
};

stRPMergeColumn *stRPMergeColumn_construct(uint64_t maskFrom, uint64_t maskTo);

void stRPMergeColumn_destruct(stRPMergeColumn *mColumn);

void stRPMergeColumn_print(stRPMergeColumn *mColumn, FILE *fileHandle, bool includeCells);

stRPMergeCell *stRPMergeColumn_getNextMergeCell(stRPCell *cell, stRPMergeColumn *mergeColumn);

stRPMergeCell *stRPMergeColumn_getPreviousMergeCell(stRPCell *cell, stRPMergeColumn *mergeColumn);

int64_t stRPMergeColumn_numberOfPartitions(stRPMergeColumn *mColumn);

/*
 * _stRPMergeCell
 * Merge cell of read partitioning hmm
 */
struct _stRPMergeCell {
	uint64_t fromPartition;
	uint64_t toPartition;
	double forwardLogProb, backwardLogProb;
};

stRPMergeCell *stRPMergeCell_construct(uint64_t fromPartition,
									   uint64_t toPartition, stRPMergeColumn *mColumn);

void stRPMergeCell_destruct(stRPMergeCell *mCell);

void stRPMergeCell_print(stRPMergeCell *mCell, FILE *fileHandle);

double stRPMergeCell_posteriorProb(stRPMergeCell *mCell, stRPMergeColumn *mColumn);

/*
 * _stGenomeFragment
 * String to represent genotype and haplotype inference from an HMM
 */
struct _stGenomeFragment {
	// The reference coordinates of the genotypes & other read info
	stReference *reference; // The reference this fragment refers to
	uint64_t refStart; // First site in the reference
	uint64_t length; // The number of sites

	// Reads
	stSet *reads1; // The reads in the first partition
	stSet *reads2; // The reads in the second partition

	// A string where each element represents the predicted genotype at the corresponding
	// position.
	// A genotype is represented by an integer in the range [0, allele_number**2), where
	// where allele_number is the number alleles at the given site
	// A genotype expresses two characters. For two characters x, y represented by two integers
	// in [0, allele_number) then the genotype is expressed as x * allele_number + y if x <= y
	// else y * allele_number + x
	uint64_t *genotypeString;

	// Strings representing the predicted haplotypes, where each element is a reference to an allele
	uint64_t *haplotypeString1;
	uint64_t *haplotypeString2;
	uint64_t *ancestorString; // predicted ancestral alleles

	// Number of reads supporting each allele
	uint64_t *readsSupportingHaplotype1;
	uint64_t *readsSupportingHaplotype2;

	// An array of genotype posterior probabilities,
	// each between 0 and 1, for the corresponding genotypes
	// in the genotype string
	float *genotypeProbs;
	float *haplotypeProbs1;
	float *haplotypeProbs2;
};

stGenomeFragment *
stGenomeFragment_constructEmpty(stReference *ref, uint64_t refStart, uint64_t length, stSet *reads1, stSet *reads2);

stGenomeFragment *stGenomeFragment_construct(stRPHmm *hmm, stList *path);

void stGenomeFragment_destruct(stGenomeFragment *genomeFragment);

void stGenomeFragment_refineGenomeFragment(stGenomeFragment *gF,
										   stRPHmm *hmm, stList *path, int64_t maxIterations);

void stGenomeFragment_printPartitionAsCSV(stGenomeFragment *gF, FILE *fh, stRPHmmParameters *params, bool hap1,
        stSet *printedReads);

void stGenomeFragment_phaseBamChunkReads(stGenomeFragment *gf, stHash *readsToPSeqs, stList *reads,
										 stSet **readsBelongingToHap1, stSet **readsBelongingToHap2,
										 stRPHmmParameters *params);

double getLogProbOfReadGivenHaplotype(const uint64_t *haplotypeString, int64_t start, int64_t length,
									  stProfileSeq *profileSeq, stReference *ref);




/*
 * Polish functions
 */

/*
 * Parameter object for polish algorithm
 */

typedef enum {
    SCM_RANDOM=0,
    SCM_SIZE_DESC=1,
} ShuffleChunksMethod;

struct _polishParams {
	bool useRunLengthEncoding;

	// Models for comparing sequences
	Alphabet *alphabet; // The alphabet object
	StateMachine *stateMachineForGenomeComparison; // Statemachine for comparing two haplotypes
	StateMachine *stateMachineForForwardStrandRead; // Statemachine for a forward strand read aligned to a reference assembly.
	StateMachine *stateMachineForReverseStrandRead; // Statemachine for a reverse strand read aligned to a reference assembly.
	PairwiseAlignmentParameters *p; // Parameters object used for aligning
	RepeatSubMatrix *repeatSubMatrix; // Repeat counts model used for predicting repeat counts of RLE sequences
	bool useRepeatCountsInAlignment; // Use repeat counts in comparing reads to a reference
	bool useReadAlleles; // Use read substrings rather than substrings (alleles) sampled from paths in the POA in the polish algorithm
	bool useReadAllelesInPhasing; // Use read substrings rather than substrings (alleles) sampled from paths in the POA in the phasing algorithm

	// chunking configuration
	bool shuffleChunks;
    ShuffleChunksMethod shuffleChunksMethod;
	bool includeSoftClipping;
	uint64_t chunkSize;
	uint64_t chunkBoundary;
	// input reads configuration
	uint64_t maxDepth;
	uint64_t excessiveDepthThreshold; // depth threshold where we randomly discard reads on initial reading
	bool includeSecondaryAlignments;
	bool includeSupplementaryAlignments;
    uint64_t filterAlignmentsWithMapQBelowThisThreshold;
    // other configuration
    double candidateVariantWeight; // The fraction (from 0 to 1) of the average position coverage needed to define a candidate variant
    uint64_t columnAnchorTrim; // The min distance between a column anchor and a candidate variant
    uint64_t maxConsensusStrings; // The maximum number of different consensus strings to consider for a substring.
    uint64_t maxPoaConsensusIterations; // Maximum number of poa_consensus / realignment iterations
    uint64_t minPoaConsensusIterations; // Minimum number of poa_consensus / realignment iterations
    uint64_t maxRealignmentPolishIterations; // Maximum number of poa_polish iterations
    uint64_t minRealignmentPolishIterations; // Minimum number of poa_polish iterations
    uint64_t filterReadsWhileHaveAtLeastThisCoverage; // Only filter read substrings if we have at least this coverage
    bool skipHaploidPolishingIfDiploid; // If doing diploid polishing, flag specifies if to skip initial diploid polish
    // at a locus
    double minAvgBaseQuality; // Minimum average base quality to include a substring for consensus finding
    double hetSubstitutionProbability; // The probability of a heterozygous variant
    double hetRunLengthSubstitutionProbability; // The probability of a heterozygous run length

    // Poa parameters
    bool poaConstructCompareRepeatCounts; // use the repeat counts in deciding if an indel can be shifted
    double referenceBasePenalty; // used by poa_getConsensus to weight against picking the reference base
    double *minPosteriorProbForAlignmentAnchors; // used by by poa_getAnchorAlignments to determine which alignment pairs
    // to use for alignment anchors during poa_realignIterative, of the form of even-length array of form
	// [ min_posterio_anchor_prob_1, diagonal_expansion_1,  min_posterio_anchor_prob_2, diagonal_expansion_2, ... ]
	int64_t minPosteriorProbForAlignmentAnchorsLength;  // Length of array minPosteriorProbForAlignmentAnchors
};

PolishParams *polishParams_readParams(FILE *fileHandle);

PolishParams *polishParams_constructEmpty();

void polishParams_printParameters(PolishParams *polishParams, FILE *fh);

void polishParams_destruct(PolishParams *polishParams);

/*
 * Basic data structures for representing a POA alignment.
 */

struct _Poa {
	Alphabet *alphabet; // The alphabet of bases
	uint64_t maxRepeatCount; // The maximum repeat count, exclusive
	RleString *refString; // The reference string, encoded using RLE
	stList *nodes;
};

struct _poaNode {
	stList *inserts; // Inserts that happen immediately after this position
	stList *deletes; // Deletes that happen immediately after this position
	char base; // Char representing base, e.g. 'A', 'C', etc.
	uint64_t repeatCount; // Repeat count of base
	double *baseWeights; // Weight given to each possible base
	double *repeatCountWeights; // Weight given to each possible repeat count
	stList *observations; // Individual events representing event, a list of PoaObservations
};

struct _poaInsert {
	RleString *insert; // RLE string representing characters of insert e.g. "GAT" with repeat counts "121", etc.
	double weightForwardStrand;
	double weightReverseStrand;
	stList *observations; // Individual events representing event, a list of PoaObservations
};

struct _poaDelete {
	int64_t length; // Length of delete
	double weightForwardStrand;
	double weightReverseStrand;
	stList *observations; // Individual events representing event, a list of PoaObservations
};

struct _poaBaseObservation {
	int64_t readNo;
	int64_t offset;
	double weight;
};


/*
 * Poa functions.
 */

double poaInsert_getWeight(PoaInsert *toInsert);

double poaDelete_getWeight(PoaDelete *toDelete);

/*
 * Creates a POA representing the given RLE reference sequence, with one node for each reference base and a
 * prefix 'N'/1 base to represent place to add inserts/deletes that precede the first position of the reference.
 */
Poa *poa_getReferenceGraph(RleString *reference, Alphabet *alphabet, uint64_t maxRepeatCount);

/*
 * Adds to given POA the matches, inserts and deletes from the alignment of the given read to the reference.
 * Adds the inserts and deletes so that they are left aligned.
 */
void poa_augment(Poa *poa, RleString *read, bool readStrand, int64_t readNo, stList *matches, stList *inserts,
				 stList *deletes,
				 PolishParams *polishParams);

/*
 * Creates a POA representing the reference and the expected inserts / deletes and substitutions from the
 * alignment of the given set of reads aligned to the reference. Anchor alignments is a set of pairwise
 * alignments between the reads and the reference sequence. There is one alignment for each read. See
 * poa_getAnchorAlignments. The anchorAlignments can be null, in which case no anchors are used.
 */
Poa *poa_realign(stList *bamChunkReads, stList *alignments, RleString *reference, PolishParams *polishParams);

/*
 * Creates a POA representing the reference and the inserts / deletes and substitutions only in the anchor
 * aligments.
 */
Poa *poa_realignOnlyAnchorAlignments(stList *bamChunkReads, stList *anchorAlignments, RleString *reference,
									 PolishParams *polishParams);

/*
 * Generates a set of anchor alignments for the reads aligned to a consensus sequence derived from the poa.
 * These anchors can be used to restrict subsequent alignments to the consensus to generate a new poa.
 * PoaToConsensusMap is a map from the positions in the poa reference sequence to the derived consensus
 * sequence. See poa_getConsensus for description of poaToConsensusMap. If poaToConsensusMap is NULL then
 * the alignment is just the reference sequence of the poa.
 */
stList *poa_getAnchorAlignments(Poa *poa, const int64_t *poaToConsensusMap, int64_t noOfReads,
								PolishParams *polishParams);

/*
 * Generates a set of maximal expected alignments for the reads aligned to the the POA reference sequence.
 * Unlike the draft anchor alignments, these are designed to be complete, high quality alignments.
 */
stList *poa_getReadAlignmentsToConsensus(Poa *poa, stList *bamChunkReads, PolishParams *polishParams);


/*
 * Prints representation of the POA.
 */
void poa_print(Poa *poa, FILE *fH,
			   stList *bamChunkReads,
			   float indelSignificanceThreshold);

/*
 * Prints a tab separated version of the POA graph.
 */
void poa_printDOT(Poa *poa, FILE *fH, stList *bamChunkReads);

/*
 * Prints a comma separated version of the POA graph.
 *
 * Format is csv with one line per reference position.
 * The CSV has the following fields (in order):
 *
 * REF_INDEX: Index of base in POA backbone, 0 is the position before the first actual base, used to represent prefix indels
 * REF_BASE: Base at REF_INDEX in the POA backbone
 * TOTAL_WEIGHT: (float) The expected read coverage of the base at REF_INDEX (float >= 0)
 * FRACTION_POS_STRAND: (float) Fraction of weight from positive strand reads.
 * for each base, c, in alphabet (bases A, C, G, T for DNA, in order):
 * FRACTION_BASE_c_WEIGHT: (float between 0 and 1) The read weight of reads that have a c aligned to the reference base, as a fraction of total weight.
 * FRACTION_BASE_c_POS_STRAND: (float between 0 and 1) The read weight of positive strand reads that have a c aligned to the reference base,
 * as a fraction of total weight of reads have a c aligned to the reference base.
 *
 * for each repeat count (from 1 to maxRepeatCount):
 * LOG_PROB_REPEAT_COUNT_i: Log probability of repeat count i (according to Bayesian model).
 *
 * for each insert:
 * INSERT_SEQ: The full insertion sequence, not in RLE space
 * TOTAL_WEIGHT: (float) The expected read coverage of the insert
 * FRACTION_POS_STRAND: (float between 0 and 1) The fraction of read weight from positive strand reads.
 *
 * for each delete:
 * DELETE_LENGTH: Length of the delete
 * TOTAL_WEIGHT: (float) The expected read coverage of the delete
 * FRACTION_POS_STRAND: (float between 0 and 1) The fraction of read weight from positive strand reads.
 */
void poa_printCSV(Poa *poa, FILE *fH,
				  stList *bamChunkReads, RepeatSubMatrix *repeatSubMatrix,
				  float indelSignificanceThreshold);

/*
 * Similar to poa_printCSV, but for a phased poa, giving weights for the two haplotypes separately.
 *
 * Format is csv with one line per reference position.
 * The CSV has the following fields (in order):
 * REF_INDEX: Index of base in POA backbone, 0 is the position before the first actual base, used to represent prefix indels
 * REF_BASE: Base at REF_INDEX in the POA backbone
 * TOTAL_WEIGHT: (float) The expected read coverage of the base at REF_INDEX (float >= 0)
 * FRACTION_HAP1_WEIGHT: The fraction (float between 0 and 1) of total weight that is in reads that are part of the first haplotype partition.
 * FRACTION_HAP2_WEIGHT: The fraction (float between 0 and 1) of total weight that is in reads that are part of the second haplotype partition.
 * FRACTION_POS_STRAND_HAP1: (float between 0 and 1) The read weight of first haplotype, positive strand reads as a fraction of all first haplotype read weight.
 * FRACTION_POS_STRAND_HAP2: (float between 0 and 1) The read weight of second haplotype, positive strand reads as a fraction of all second haplotype read weight.
 * for each base, c, in alphabet (bases A, C, G, T for DNA, in order):
 * NORM_BASE_c_WEIGHT: (float between 0 and 1) The read weight of reads that have a c aligned to the reference base, as a fraction of total weight.
 * FRACTION_BASE_c_HAP1: (float between 0 and 1) The read weight of first haplotype reads that have a c aligned to the reference base, as a fraction of total weight for first haplotype reads.
 * FRACTION_BASE_c_HAP2 (float between 0 and 1) The read weight of second haplotype reads that have a c aligned to the reference base, as a fraction of total weight for second haplotype reads.
 * FRACTION_BASE_c_POS_STRAND_HAP1: (float between 0 and 1) The read weight of positive-strand first haplotype reads that have a c aligned to the reference base,
 * as a fraction of total weight for first haplotype reads with a c aligned to the reference base.
 * FRACTION_BASE_c_POS_STRAND_HAP2: (float between 0 and 1) The read weight of positive-strand second haplotype reads that have a c aligned to the reference base,
 * as a fraction of total weight for second haplotype reads with a c aligned to the reference base.
 *
 * for each repeat count for first haplotype (from 1 to maxRepeatCount):
 * LOG_PROB_HAP1_REPEAT_COUNT_i: Log probability of repeat count i for haplotype1 (according to Bayesian model).
 *
 * for each repeat count for first haplotype (from 1 to maxRepeatCount):
 * LOG_PROB_HAP2_REPEAT_COUNT_i: Log probability of repeat count i for haplotype2 (according to Bayesian model).
 *
 * for each insert:
 * INSERT_SEQ: The full insertion sequence, not in RLE space
 * TOTAL_WEIGHT: (float) The expected read coverage of the insert
 * FRACTION_HAP1_WEIGHT: The fraction (float between 0 and 1) of total weight that is in reads that are part of the first haplotype partition.
 * FRACTION_HAP2_WEIGHT: The fraction (float between 0 and 1) of total weight that is in reads that are part of the second haplotype partition.
 * FRACTION_POS_STRAND_HAP1: (float between 0 and 1) The read weight of first haplotype, positive strand reads as a fraction of all first haplotype read weight.
 * FRACTION_POS_STRAND_HAP2: (float between 0 and 1) The read weight of second haplotype, positive strand reads as a fraction of all second haplotype read weight.
 *
 * for each delete:
 * DELETE_LENGTH: Length of the insert:
 * TOTAL_WEIGHT: (float) The expected read coverage of the insert
 * FRACTION_HAP1_WEIGHT: The fraction (float between 0 and 1) of total weight that is in reads that are part of the first haplotype partition.
 * FRACTION_HAP2_WEIGHT: The fraction (float between 0 and 1) of total weight that is in reads that are part of the second haplotype partition.
 * FRACTION_POS_STRAND_HAP1: (float between 0 and 1) The read weight of first haplotype, positive strand reads as a fraction of all first haplotype read weight.
 * FRACTION_POS_STRAND_HAP2: (float between 0 and 1) The read weight of second haplotype, positive strand reads as a fraction of all second haplotype read weight.
 */
void poa_printPhasedCSV(Poa *poa, FILE *fH,
						stList *bamChunkReads, stSet *readsInHap1, stSet *readsInHap2,
						RepeatSubMatrix *repeatSubMatrix,
						float indelSignificanceThreshold);

/*
 * Print individual repeat count observations.
 */
void poa_printRepeatCountsCSV(Poa *poa, FILE *fH, stList *bamChunkReads);

/*
 * Prints some summary stats on the POA.
 */
void poa_printSummaryStats(Poa *poa, FILE *fH);

/*
 * Creates a consensus reference sequence from the POA. poaToConsensusMap is a pointer to an
 * array of integers of length poa->refString->length, giving the index of the reference positions
 * alignment to the consensus sequence, or -1 if not aligned. It is initialised as a
 * return value of the function.
 */
RleString *poa_getConsensus(Poa *poa, int64_t **poaToConsensusMap, PolishParams *polishParams);

RleString *poa_polish(Poa *poa, stList *bamChunkReads, PolishParams *params,
					  int64_t **poaToConsensusMap);


/*
 * Iteratively uses poa_getConsensus and poa_polish to refine the median reference sequence
 * for the given reads and the starting reference.
 *
 * Allows the specification of the min and max number of realignment cycles.
 */
Poa *poa_realignIterative(Poa *poa, stList *bamChunkReads,
						  PolishParams *polishParams, bool hmmNotRealign,
						  int64_t minIterations, int64_t maxIterations);

/*
 * Convenience function that iteratively polishes sequence using poa_getConsensus and then poa_polish for
 * a specified number of iterations.
 */
Poa *poa_realignAll(stList *bamChunkReads, stList *anchorAlignments, RleString *reference,
					PolishParams *polishParams);

/*
 * Greedily evaluate the top scoring indels.
 */
Poa *poa_checkMajorIndelEditsGreedily(Poa *poa, stList *bamChunkReads, PolishParams *polishParams);

void poa_destruct(Poa *poa);

/*
 * Compare functions
 */
int cmpAlignedPairsByCoordinates(const void *a, const void *b);
int cmpAlignedPairsByInvertedCoordinates(const void *a, const void *b);


double *poaNode_getStrandSpecificBaseWeights(PoaNode *node, stList *bamChunkReads,
											 double *totalWeight, double *totalPositiveWeight,
											 double *totalNegativeWeight, Alphabet *a, stSet *readsToInclude);

/*
 * Finds shift, expressed as a reference coordinate, that the given substring str can
 * be shifted left in the refString, starting from a match at refStart.
 *
 * If compareRepeatCounts is non-zero then comparison includes exact comparison of repeat counts.
 */
int64_t getShift(RleString *refString, int64_t refStart, RleString *str, bool compareRepeatCounts);

/*
 * Get sum of weights for reference bases in poa - proxy to agreement of reads
 * with reference.
 */
double poa_getReferenceNodeTotalMatchWeight(Poa *poa);

/*
 * Get sum of weights for delete in poa - proxy to delete disagreement of reads
 * with reference.
 */
double poa_getDeleteTotalWeight(Poa *poa);

/*
 * Get sum of weights for inserts in poa - proxy to insert disagreement of reads
 * with reference.
 */
double poa_getInsertTotalWeight(Poa *poa);

/*
 * Get sum of weights for non-reference bases in poa - proxy to disagreement of read positions
 * aligned with reference.
 */
double poa_getReferenceNodeTotalDisagreementWeight(Poa *poa);

/*
 * Reestimates the repeat counts of the poa backbone using the Bayesian model, repeatSubMatrix.
 * Changes the repeat counts of the backbone bases in place
 */
void poa_estimateRepeatCountsUsingBayesianModel(Poa *poa, stList *bamChunkReads, RepeatSubMatrix *repeatSubMatrix);

/*
 * As poa_estimateRepeatCountsUsingBayesianModel, but using a phasing.
 */
void poa_estimatePhasedRepeatCountsUsingBayesianModel(Poa *poa, stList *bamChunkReads,
													  RepeatSubMatrix *repeatSubMatrix, stSet *readsBelongingToHap1,
													  stSet *readsBelongingToHap2, PolishParams *params);

/*
 * Uses phasing to estimate ML bases.
 */
void poa_estimatePhasedBasesUsingBayesianModel(Poa *poa, stList *bamChunkReads, stSet *readsBelongingToHap1,
											   PolishParams *params);

// Data structure for representing RLE strings
struct _rleString {
	char *rleString; //Run-length-encoded (RLE) string
	uint64_t *repeatCounts; // Count of repeat for each position in rleString
	uint64_t length; // Length of the rleString
	uint64_t nonRleLength; // Length of the expanded, non-rle string
};

/*
 * Returns a string "cXrepeatCount", e.g. c='a', repeatCount=4 returns "aaaa".
 */
char *expandChar(char c, uint64_t repeatCount);

RleString *rleString_construct(char *string);

RleString *rleString_constructPreComputed(char *rleChars, const uint8_t *rleCounts);

RleString *rleString_construct_no_rle(char *string);

void rleString_destruct(RleString *rlString);

RleString *rleString_copy(RleString *rleString);

RleString *rleString_copySubstring(RleString *rleString, uint64_t start, uint64_t length);

bool rleString_eq(RleString *r1, RleString *r2);

/*
 * Debug output friendly version of rleString on one line.
 * Does not print any new lines.
 */
void rleString_print(RleString *rleString, FILE *f);

/*
 * Cyclic rotatation of the rle string so that the suffix of str of rotationLength is removed and made the prefix
 * of the string.
 */
void rleString_rotateString(RleString *str, int64_t rotationLength, bool mergeEnds);

/*
 * Generates the expanded non-rle string.
 */
char *rleString_expand(RleString *rleString);

/*
 * Gets a symbol sub-string from a given RLE string.
 */
SymbolString rleString_constructSymbolString(RleString *s, int64_t start, int64_t length, Alphabet *a, bool includeRepeatCounts, uint64_t maxRunLengthExclusive);

/*
 * Gets an array giving the position in the rleString of a corresponding position in the expanded string.
 */
uint64_t *rleString_getNonRleToRleCoordinateMap(RleString *rleString);

uint64_t *rleString_getRleToNonRleCoordinateMap(RleString *rleString);

uint8_t *rleString_rleQualities(RleString *rleString, const uint8_t *qualities);

// Data structure for storing log-probabilities of observing
// one repeat count given another
struct _repeatSubMatrix {
	Alphabet *alphabet;
	double *repeatCountSubstitutionLogProb;
	double *baseLogProbs_AT;
	double *baseLogProbs_GC;
	double *logProbabilities;
	int64_t maximumRepeatLength; // The maximum repeat count length, exclusive
	int64_t maxEntry;
};

int64_t getMax(double *values, int64_t length,
			   double *maxValue);

RepeatSubMatrix *repeatSubMatrix_constructEmpty(Alphabet *alphabet);

void repeatSubMatrix_destruct(RepeatSubMatrix *repeatSubMatrix);

/*
 * Gets the log probability of observing a given repeat conditioned on an underlying repeat count and base.
 */
double repeatSubMatrix_getLogProb(RepeatSubMatrix *repeatSubMatrix, Symbol base, bool strand,
								  int64_t observedRepeatCount, int64_t underlyingRepeatCount);

/*
 * As gets, but returns the address.
 */
double *
repeatSubMatrix_setLogProb(RepeatSubMatrix *repeatSubMatrix, Symbol base, bool strand, int64_t observedRepeatCount,
						   int64_t underlyingRepeatCount);

/*
 * Gets the log probability of observing a given set of repeat observations conditioned on an underlying repeat count and base.
 */
double repeatSubMatrix_getLogProbForGivenRepeatCount(RepeatSubMatrix *repeatSubMatrix, Symbol base,
													 stList *observations, stList *bamChunkReads,
													 int64_t underlyingRepeatCount);

/*
 * Gets the maximum likelihood underlying repeat count for a given set of observed read repeat counts.
 * Puts the ml log probility in *logProbabilty.
 */
int64_t repeatSubMatrix_getMLRepeatCount(RepeatSubMatrix *repeatSubMatrix, Symbol base, stList *observations,
										 stList *bamChunkReads, double *logProbability);

/*
 * Get the log probabilities of repeat counts from minRepeatLength (inclusive) to maxRepeatLength (inclusive)
 */
void repeatSubMatrix_getRepeatCountProbs(RepeatSubMatrix *repeatSubMatrix, Symbol base, stList *observations,
										 stList *bamChunkReads, double *logProbabilities, int64_t minRepeatLength,
										 int64_t maxRepeatLength);

/*
 * As repeatSubMatrix_getMLRepeatCount, but for a phasing of the reads.
 */
int64_t
repeatSubMatrix_getPhasedMLRepeatCount(RepeatSubMatrix *repeatSubMatrix, int64_t existingRepeatCount, Symbol base,
									   stList *observations,
									   stList *bamChunkReads, double *logProbability, stSet *readsBelongingToHap1,
									   stSet *readsBelongingToHap2, PolishParams *params);

/*
 * Get the minimum and maximum repeat count observations.
 */
void repeatSubMatrix_getMinAndMaxRepeatCountObservations(RepeatSubMatrix *repeatSubMatrix, stList *observations,
														 stList *bamChunkReads, int64_t *minRepeatLength,
														 int64_t *maxRepeatLength);

/*
 * Translate a sequence of aligned pairs (as stIntTuples) whose coordinates are monotonically increasing
 * in both underlying sequences (seqX and seqY) into an equivalent run-length encoded space alignment.
 */
stList *runLengthEncodeAlignment(stList *alignment,
								 const uint64_t *seqXNonRleToRleCoordinateMap,
								 const uint64_t *seqYNonRleToRleCoordinateMap);

/*
 * Make edited string with given insert. Edit start is the index of the position to insert the string.
 */
char *addInsert(char *string, char *insertString, int64_t editStart);

/*
 * Make edited string with given insert. Edit start is the index of the first position to delete from the string.
 */
char *removeDelete(char *string, int64_t deleteLength, int64_t editStart);

/*
 * Generates aligned pairs and indel probs, but first crops reference to only include sequence from first
 * to last anchor position.
 */
void getAlignedPairsWithIndelsCroppingReference(RleString *reference,
												RleString *read, bool readStrand, stList *anchorPairs,
												stList **matches, stList **inserts, stList **deletes,
												PolishParams *polishParams);

/*
 * Functions for processing BAMs
 */

// TODO: MOVE BAMCHUNKER TO PARSER .c

typedef struct _bamChunker {
	// file locations
	char *bamFile;
	// configuration
	uint64_t chunkSize;
	uint64_t chunkBoundary;
	bool includeSoftClip;
	PolishParams *params;
	// internal data
	stList *chunks;
	uint64_t chunkCount;
	stHash *readEnumerator;
} BamChunker;

typedef struct _bamChunk {
    char *refSeqName;          // name of contig
    int64_t chunkIdx;
    int64_t chunkOverlapStart;  // the first 'position' where we have an aligned read
    int64_t chunkStart;        // the actual boundary of the chunk, calculations from chunkMarginStart to chunkStart
    //  should be used to initialize the probabilities at chunkStart
    int64_t chunkEnd;          // same for chunk end
    int64_t chunkOverlapEnd;    // no reads should start after this position
    int64_t estimatedDepth;		// for ranking chunk order
    BamChunker *parent;        // reference to parent (may not be needed)
} BamChunk;

typedef struct _bamChunkSubstrings {
    stList *readSubstrings;
    stList *readSubstringQualities;
    stList *vcfEntries;
} BamChunkReadVcfEntrySubstrings;

typedef struct _bamChunkRead {
	char *readName;            // read name
	RleString *rleRead;        // rle read
	uint8_t *qualities;            // quality scores. will be NULL if not given, else will be of length rleRead->length
	bool forwardStrand;            // whether the alignment is matched to the forward strand
	int64_t fullReadLength;   // total length for whole read (not just chunk portion)
	BamChunkReadVcfEntrySubstrings *bamChunkReadVcfEntrySubstrings; // for ultra-fast phasing work
} BamChunkRead;

BamChunkRead *bamChunkRead_construct2(char *readName, char *nucleotides, uint8_t *qualities, bool forwardStrand,
                                      bool useRunLengthEncoding);
BamChunkRead *bamChunkRead_construct3(char *readName, char *nucleotides, uint8_t *qualities, bool forwardStrand,
                                      int64_t fullReadLength, bool useRunLengthEncoding);
BamChunkRead *bamChunkRead_constructWithVcfEntrySubstrings(char *readName, bool forwardStrand, int64_t fullReadLength,
        BamChunkReadVcfEntrySubstrings *bcrves);

BamChunkRead *bamChunkRead_constructCopy(BamChunkRead *copy);

void bamChunkRead_destruct(BamChunkRead *bamChunkRead);

BamChunkReadVcfEntrySubstrings *bamChunkReadVcfEntrySubstrings_construct();
BamChunkReadVcfEntrySubstrings *bamChunkReadVcfEntrySubstrings_construct2( stList *readSubstrings,
                                                                           stList *readSubstringQualities,
                                                                           stList *vcfEntries);
BamChunkReadVcfEntrySubstrings *bamChunkReadVcfEntrySubstrings_constructCopy(BamChunkReadVcfEntrySubstrings *bcrves);

void bamChunkReadVcfEntrySubstrings_destruct(BamChunkReadVcfEntrySubstrings *bcrs);

void bamChunkReadVcfEntrySubstrings_saveSubstring(BamChunkReadVcfEntrySubstrings *bcrs, char *read, uint8_t *qualities,
                                                  VcfEntry *vcfEntry);


/*
 * Generates the expanded non-rle version of bam chunk read nucleotide sequence.
 */
char *bamChunkRead_rleExpand(BamChunkRead *read);

typedef struct _bamChunkReadSubstring {
	BamChunkRead *read; // The parent read from which the substring arises from
	int64_t start; // The 0 based offset of the start position in the parent read (inclusive)
	int64_t length; // The length of the substring
	double qualValue;
	RleString *substring; // for the ultra-fast version where we don't use the BamChunkRead to get the substring
} BamChunkReadSubstring;

/*
 * Gets the RLE substring for the bam chunk read substring.
 */
RleString *bamChunkReadSubstring_getRleString(BamChunkReadSubstring *readSubstring);

/*
 * Remove overlap between two overlapping strings. Returns max weight of split point.
 */
int64_t removeOverlap(char *prefixString, int64_t prefixStringLength, char *suffixString, int64_t suffixStringLength,
					  int64_t approxOverlap, PolishParams *polishParams,
					  int64_t *prefixStringCropEnd, int64_t *suffixStringCropStart);

/*
 * View functions
 */

struct _refMsaView {
	int64_t refLength; // The length of the reference sequence
	char *refSeq; // The reference sequence - this is not copied by the constructor
	char *refSeqName; // The reference sequence name - this is not copied by the constructor, and can be NULL
	int64_t seqNo; // The number of non-ref sequences aligned to the reference
	stList *seqs; // The non-ref sequences - - this is not copied by the constructor
	stList *seqNames; // The non-ref sequence names - this is not copied by the constructor, and can be NULL
	int64_t *seqCoordinates; // A matrix giving the coordinates of the non-reference sequence
	// as aligned to the reference sequence
	int64_t *maxPrecedingInsertLengths; // The maximum length of an insert in
	// any of the sequences preceding the reference positions
	int64_t **precedingInsertCoverages; // The number of sequences with each given indel position
};

/*
 * Get the coordinate in the given sequence aligned to the given reference position. Returns -1 if
 * no sequence position is aligned to the reference position.
 */
int64_t msaView_getSeqCoordinate(MsaView *view, int64_t refCoordinate, int64_t seqIndex);

/*
 * Gets the length of any insert in the given sequence preceding the given reference position. If
 * the sequence is not aligned at the given reference position returns 0.
 */
int64_t msaView_getPrecedingInsertLength(MsaView *view, int64_t rightRefCoordinate, int64_t seqIndex);

/*
 * Gets the first position in the sequence of an insert preceding the given reference position. If there
 * is no such insert returns -1.
 */
int64_t msaView_getPrecedingInsertStart(MsaView *view, int64_t rightRefCoordinate, int64_t seqIndex);

/*
 * Gets the maximum length of an indel preceding the given reference position
 */
int64_t msaView_getMaxPrecedingInsertLength(MsaView *view, int64_t rightRefCoordinate);

/*
 * Get the number of sequences with an insertion at a given position. IndelOffset if the position, from 0, of the
 indel from left-to-right.
 */
int64_t msaView_getPrecedingCoverageDepth(MsaView *view, int64_t rightRefCoordinate, int64_t indelOffset);

/*
 * Get the maximum length of an insertion at a given position with a minimum of reads supporting it.
 */
int64_t
msaView_getMaxPrecedingInsertLengthWithGivenCoverage(MsaView *view, int64_t rightRefCoordinate, int64_t minCoverage);

/*
 * Builds an MSA view for the given reference and aligned sequences.
 * Does not copy the strings or string names, just holds references.
 */
MsaView *msaView_construct(char *refSeq, char *refName,
						   stList *refToSeqAlignments, stList *seqs, stList *seqNames);

void msaView_destruct(MsaView *view);

/*
 * Prints a quick view of the MSA for debugging/browsing.
 */
void msaView_print(MsaView *view, int64_t minInsertCoverage, FILE *fh);

/*
 * Prints the repeat counts of the MSA.
 */
void msaView_printRepeatCounts(MsaView *view, int64_t minInsertCoverage,
							   RleString *refString, stList *rleStrings, FILE *fh);

/*
 * Bubble graphs
 */

typedef struct _bubble {
	uint64_t refStart; //First inclusive position
	uint64_t bubbleLength; //length of reference between bubble anchors
	RleString *refAllele; // The current reference allele
	uint64_t alleleNo; // Number of alleles
	RleString **alleles; // Array of allele strings, each an RLE string
	uint64_t readNo; // Number of reads overlapping bubble
	BamChunkReadSubstring **reads; // Array of read substrings aligned to the bubble
	float *alleleReadSupports; // An array of log-likelihoods giving the support of
	// each allele for each read, stored as [i * readNo + j], where i is the allele index
	// and j is the index of the read
	uint64_t alleleOffset; // The index of the first allele in this bubble
	// in a sequence of all alleles in the bubble graph, ordered first by bubble then
	// by order in the bubble.
	stList *variantPositionOffsets;
} Bubble;

typedef struct _bubbleGraph {
	RleString *refString; // The reference string
	uint64_t bubbleNo; // The number of bubbles
	Bubble *bubbles; // An array of bubbles
	uint64_t totalAlleles; // Sum of alleles across bubbles
} BubbleGraph;

/*
 * Get a consensus path through bubble graph by picking the highest
 * likelihood allele at each bubble. Returned as a string of bg->refLength integers,
 * each denoting a chosen allele.
 */
uint64_t *bubbleGraph_getConsensusPath(BubbleGraph *bg, PolishParams *polishParams);

/*
 * Get a consensus sequences from the bubble graph by picking the highest
 * likelihood allele at each bubble.
 */
RleString *bubbleGraph_getConsensusString(BubbleGraph *bg, uint64_t *consensusPath, int64_t **poaToConsensusMap,
										  PolishParams *polishParams);

/*
 * Create a bubble graph from a POA.
 */
BubbleGraph *bubbleGraph_constructFromPoa(Poa *poa, stList *bamChunkReads, PolishParams *params);
BubbleGraph *bubbleGraph_constructFromPoa2(Poa *poa, stList *bamChunkReads, PolishParams *params, bool phasing);

void bubbleGraph_destruct(BubbleGraph *bg);

/*
 * The the index in b->alleles of the allele with highest likelihood
 * given the reads
 */
int64_t bubble_getIndexOfHighestLikelihoodAllele(Bubble *b, PolishParams *p);

/*
 * Gets the likelihood of a given allele giving rise to the reads.
 */
double bubble_getLogLikelihoodOfAllele(Bubble *b, int64_t allele, PolishParams *p);

/*
 * Gets a set of profile sequences for the reads aligned to the bubble graph.
 * Allows them to then be phased. Returns as hash of bamChunkReads to profileSeqs.
 */
stHash *bubbleGraph_getProfileSeqs(BubbleGraph *bg, stReference *ref);

/*
 * Gets an stReference that can be used for phasing.
 */
stReference *bubbleGraph_getReference(BubbleGraph *bg, char *refName, Params *params);

/*
 * Phase bubble graph.
 */
stGenomeFragment *
bubbleGraph_phaseBubbleGraph(BubbleGraph *bg, stReference *ref, stList *reads, Params *params, stHash **readsToPSeqs);

/*
 * Get Poa from bubble graph.
 */
Poa *bubbleGraph_getNewPoa(BubbleGraph *bg, uint64_t *consensusPath, Poa *poa, stList *reads, Params *params);

/*
 * Gets the strand support skew for each allele.
 */
void bubble_calculateStrandSkews(Bubble *b, double *skews);

/*
 * Gets the p-value for the bubble having a phased strand-skew
 */
double bubble_phasedStrandSkew(Bubble *b, stHash *readsToPSeqs, stGenomeFragment *gf);

/*
 * Returns fraction of bubbles with significant allele-strand phase skew.
 */
double bubbleGraph_skewedBubbles(BubbleGraph *bg, stHash *readsToPSeqs, stGenomeFragment *gf);

/*
 * Functions for scoring strand skews using binomial model
 */
double binomialPValue(int64_t n, int64_t k);

uint128_t bionomialCoefficient(int64_t n, int64_t k);

/*
 * Phases reads based on existing partition
 */
BubbleGraph *bubbleGraph_partitionFilteredReads(Poa *poa, stList *bamChunkReads, stGenomeFragment *gF,
												BubbleGraph *bg, BamChunk *bamChunk, uint64_t *reference_rleToNonRleCoordMap,
                                                stSet *hap1Reads, stSet *hap2Reads, PolishParams *params, FILE *out,
                                                char *logIdentifier);

/*
 * Phases reads from existing partition based on vcfEntries
 */
BubbleGraph *bubbleGraph_partitionFilteredReadsFromVcfEntries(stList *bamChunkReads, stGenomeFragment *gF,
                                                              BubbleGraph *bg, stList *vcfEntriesToBubbles, stSet *hap1Reads,
                                                              stSet *hap2Reads, Params *params, char *logIdentifier);

/*
 * For tracking Bubble Graph stuff
 */
void bubbleGraph_saveBubblePhasingInfo(BamChunk *bamChunk, BubbleGraph *bg, stHash *readsToPSeqs, stGenomeFragment *gF,
        uint64_t *reference_rleToNonRleCoordMap, FILE *out);

/*
 * For logging while multithreading
 */
char *getLogIdentifier();

/*
 * Run length encoded model emission model for state machine that uses a repeat sub matrix.
 */
Emissions *rleNucleotideEmissions_construct(Emissions *emissions, RepeatSubMatrix *repeatSubMatrix, bool strand);

/*
 * ChunkToStitch
 */

typedef struct _chunkToStitch {
    /*
     * Object for managing the output of a polished sequence.
     */
    bool startOfSequence; // Indicates if it is the first chunk in a sequence
    int64_t chunkOrdinal; // The index of the chunk in the sequence of chunks

    char *seqName; // The name of the sequence being polished

    char *seqHap1; // The primary sequence
    char *seqHap2; // The secondary haplotype, may be null.

    // Following from the output CSV files, each line corresponds to a position in the sequences
    // Each can be null.

    // Lines (strings) from the POA:
    stList *poaHap1StringsLines; // If not diploid, this is used to store the POA lines
    stList *poaHap2StringsLines;

    // Lines from the repeat count file
    stList *repeatCountLinesHap1; // If not diploid, this is used to store the POA lines
    stList *repeatCountLinesHap2;

    // Both these will be present if phasing
    stList *readsHap1Lines; // Reads from primary sequence
    stList *readsHap2Lines; // Reads from the secondary sequence

    // For tracking which chunks were stitched
    bool wasSwitched;
} ChunkToStitch;

ChunkToStitch *chunkToStitch_construct(char *seqName, int64_t chunkOrdinal, bool phased,
                                       bool initRepeatCounts, bool initPoa);
void chunkToStitch_destruct(ChunkToStitch *chunkToStitch);

/*
 * Stitching code
 */

OutputChunkers *
outputChunkers_construct(int64_t noOfOutputChunkers, Params *params, char *outputSequenceFile, char *outputPoaFile,
						 char *outputReadPartitionFile, char *outputRepeatCountFile, char *hap1Suffix, char *hap2Suffix,
						 bool inMemoryBuffers);

void
outputChunkers_processChunkSequence(OutputChunkers *outputChunkers, int64_t chunker, int64_t chunkOrdinal,
									char *sequenceName, Poa *poa,
									stList *reads);

void outputChunkers_processChunkSequencePhased(OutputChunkers *outputChunkers, int64_t chunker, int64_t chunkOrdinal,
											   char *sequenceName, Poa *poaHap1, Poa *poaHap2, stList *reads,
											   stSet *readsBelongingToHap1, stSet *readsBelongingToHap2,
											   stGenomeFragment *gF, Params *params);

void outputChunkers_stitch(OutputChunkers *outputChunkers, bool phased, int64_t chunkCount);
void outputChunkers_stitchAndTrackExtraData(OutputChunkers *outputChunkers, bool phased, int64_t chunkCount,
                                            stList *readIdsHap1, stList *readIdsHap2, bool* switchedState);

void outputChunkers_stitchLinear(OutputChunkers *outputChunkers, bool phased, Params *params);

void outputChunkers_destruct(OutputChunkers *outputChunkers);

ChunkToStitch *mergeContigChunkz(ChunkToStitch **chunks, int64_t startIdx, int64_t endIdxExclusive, bool phased,
                                 Params *params);
ChunkToStitch *mergeContigChunkzThreaded(ChunkToStitch **chunks, int64_t startIdx, int64_t endIdxExclusive, int64_t numThreads,
                                         bool phased, Params *params, char *referenceSequenceName);


/*
 * VCF data structures
 */

struct _vcfEntry {
    char *refSeqName;
    int64_t refPos;
    int64_t rawRefPosInformativeOnly;
    double quality;
    stList *alleles; //refAllele is alleles[0]
    // used by margin, include expansion around alele
    stList *alleleSubstrings;
    int64_t refAlnStart;
    int64_t refAlnStopIncl;
    VcfEntry *rootVcfEntry;
    // reports initial genotypes, updated after margin runs
    int64_t gt1;
    int64_t gt2;
    // margin calculations
    stList *alleleIdxToReads;
    double genotypeProb;
    double haplotype1Prob;
    double haplotype2Prob;
};

// vcf functions
VcfEntry *vcfEntry_construct(char *refSeqName, int64_t refPos, int64_t rawRefPos, double phredQuality,
        stList *alleles, int64_t gt1, int64_t gt2);
void vcfEntry_destruct(VcfEntry *vcfEntry);
RleString *getVcfEntryAlleleH1(VcfEntry *vcfEntry);
RleString *getVcfEntryAlleleH2(VcfEntry *vcfEntry);
stHash *parseVcf(char *vcfFile, Params *params);
stHash *parseVcf2(char *vcfFile, char *regionStr, Params *params);
int64_t binarySearchVcfListForFirstIndexAtOrAfterRefPos(stList *vcfEntries, int64_t refPos); // just exposed for testing
stList *getVcfEntriesForRegion(stHash *vcfEntries, uint64_t *rleMap, char *refSeqName, int64_t startPos,
        int64_t endPos, Params *params);
stList *getAlleleSubstrings2(VcfEntry *entry, char *referenceSeq, int64_t refSeqLen, int64_t *refStartPos,
        int64_t *refEndPosIncl, bool putRefPosInPOASpace, int64_t expansion, bool useRunLengthEncoding);
stList *getAlleleSubstrings(VcfEntry *entry, RleString *referenceSeq, Params *params,
                            int64_t *refStartPos, int64_t *refEndPosIncl, bool refPosInPOASpace);
void updateVcfEntriesWithSubstringsAndPositions(stList *vcfEntries, char *referenceSeq, int64_t refSeqLen,
        bool refPosInPOASpace, Params *params);
void updateOriginalVcfEntriesWithBubbleData(BamChunk *bamChunk, stList *bamChunkReads, stHash *readIdToIdx,
		stGenomeFragment *gF, BubbleGraph *bg, stList *chunkVcfEntriesToBubbles, stSet *hap1Reads, stSet *hap2Reads,
		char *logIdentifier);
void updateHaplotypeSwitchingInVcfEntries(BamChunker *chunker, bool *chunkWasSwitched, stHash *vcfEntryMap);
void writePhasedVcf(char *inputVcfFile, char *regionStr, char *outputVcfFile, char *phaseSetBedFile,
        stHash *vcfEntryMap, Params *params);

// bubble functions using vcfs
BubbleGraph *bubbleGraph_constructFromPoaAndVCF(Poa *poa, stList *bamChunkReads, stList *vcfEntries,
                                                PolishParams *params, bool phasing);
BubbleGraph *bubbleGraph_constructFromPoaAndVCFOnlyVCFAllele(Poa *poa, stList *bamChunkReads,
															 RleString *referenceSeq, stList *vcfEntries, Params *params);
BubbleGraph *bubbleGraph_constructFromVCFAndBamChunkReadVcfEntrySubstrings(stList *bamChunkReads, stList *vcfEntries,
                                                                           Params *params, stList **vcfEntriesToBubbleIdx);

/*
 * Misc
 */
char *getTimeDescriptorFromSeconds(int64_t seconds);

stHash *parseReferenceSequences(char *referenceFastaFile);

char *getFileBase(char *base, char *defawlt);

FILE *safe_fopen(char *file, char *openStr);

RleString *bamChunk_getReferenceSubstring(BamChunk *bamChunk, char *referenceFile, Params *params);
uint64_t *getPaddedHaplotypeString(uint64_t *hap, stGenomeFragment *gf, BubbleGraph *bg, Params *params);
stSet *bamChunkRead_to_readName(stSet *bamChunkReads);
stList *copyListOfIntTuples(stList *toCopy);
double fromLog(double theLog);
double toPhred(double prob);
double fromPhred(double phred);
void writePhasedReadInfoJSON(BamChunk *bamChunk, stList *primaryReads, stList *primaryAlignments, stList *filteredReads,
                             stList *filteredAlignments, stSet *readsInHap1, stSet *readsInHap2,
                             uint64_t *reference_rleToNonRleCoordMap, FILE *out);
void removeReadsStartingAfterChunkEnd(BamChunk *bamChunk, stList *reads, stList *alignments, char *logIdentifier);
void removeReadsOnlyInChunkBoundary(BamChunk *bamChunk, stList *reads, stList *alignments, char *logIdentifier);
stList *produceVcfEntriesFromBubbleGraph(BamChunk *bamChunk, BubbleGraph *bg, stHash *readsToPSeqs,
										 stGenomeFragment *gF, double strandSkewThreshold,
										 double readSkewThreshold);

#define CHUNK_TRUTH_READ_ID "CTRID"
#define CHUNK_TRUTH_READ_ID_LEN 5
#define CHUNK_TRUTH_READ_ID_SEP "."
typedef struct _chunkTruthHaplotypes ChunkTruthHaplotypes;
struct _chunkTruthHaplotypes {
    int64_t chunkIdx;
    stList *hap1Reads;
    stList *hap2Reads;
};
ChunkTruthHaplotypes **chunkTruthHaplotypes_construct(int64_t length);
void chunkTruthHaplotypes_destruct(ChunkTruthHaplotypes **chunkTruthHaplotypes, int64_t length);
void *chunkTruthHaplotypes_print(stList *readsInHap1, stList *readsInHap2, stList *bamChunks, int64_t length, char *filename);
void chunkTruthHaplotypes_addTruthReadsToFilteredReadSet(BamChunk *bamChunk, BamChunker *bamChunker,
        stList *readsToAdd, stList *alignmentsToAdd, RleString *rleReference, Params *params, char *logIdentifier);

/*
 * HELEN Features
 */

typedef enum {
	HFEAT_NONE=0,
	HFEAT_SIMPLE_WEIGHT=1,
	HFEAT_SPLIT_RLE_WEIGHT=2,
	HFEAT_CHANNEL_RLE_WEIGHT=3,
} HelenFeatureType;

#define POAFEATURE_SPLIT_MAX_RUN_LENGTH_DEFAULT 10
#define POAFEATURE_CHANNEL_MAX_RUN_LENGTH_DEFAULT 10

/*
 * HTS integration
 */


BamChunker *bamChunker_construct(char *bamFile, PolishParams *params);

BamChunker *bamChunker_construct2(char *bamFile, char *region, PolishParams *params, bool recordFilteredReads);

BamChunker *bamChunker_constructFromFasta(char *fastaFile, char *bamFile, char *regionStr, PolishParams *params);

BamChunker *bamChunker_copyConstruct(BamChunker *toCopy);

void bamChunker_destruct(BamChunker *bamChunker);

BamChunk *bamChunker_getChunk(BamChunker *bamChunker, int64_t chunkIdx);

BamChunk *bamChunk_construct();
BamChunk *bamChunk_construct2(char *refSeqName, int64_t chunkIndex, int64_t chunkOverlapStart, int64_t chunkStart, int64_t chunkEnd,
                              int64_t chunkOverlapEnd, int64_t depth, BamChunker *parent);

int compareBamChunkDepthByIndexInList(const void *a, const void *b, const void *chunkList);

BamChunk *bamChunk_copyConstruct(BamChunk *toCopy);

void bamChunk_destruct(BamChunk *bamChunk);


/*
 * Converts chunk of aligned reads into list of reads and alignments.
 */
uint32_t convertToReadsAndAlignments(BamChunk *bamChunk, RleString *reference, stList *reads, stList *alignments,
                                     PolishParams *polishParams);
uint32_t convertToReadsAndAlignmentsWithFiltered(BamChunk *bamChunk, RleString *reference, stList *reads,
                                                 stList *alignments, stList *filteredReads, stList *filteredAlignments,
                                                 PolishParams *polishParams);
uint32_t extractReadSubstringsAtVariantPositions(BamChunk *bamChunk, stList *vcfEntries, stList *reads,
                                                 stList *filteredReads, PolishParams *polishParams);

bool downsampleViaReadLikelihood(int64_t intendedDepth, BamChunk *bamChunk, stList *inputReads, stList *inputAlignments,
                                 stList *maintainedReads, stList *maintainedAlignments, stList *discardedReads,
                                 stList *discardedAlignments);
bool downsampleViaHetSpanLikelihood(int64_t intendedDepth, BamChunk *bamChunk, stList *vcfEntries,
                                    stList *inputReads, stList *inputAlignments, stList *maintainedReads,
                                    stList *maintainedAlignments, stList *discardedReads, stList *discardedAlignments);
bool downsampleViaFullReadLengthLikelihood(int64_t intendedDepth, BamChunk *bamChunk, stList *inputReads,
                                           stList *inputAlignments, stList *maintainedReads, stList *maintainedAlignments,
                                           stList *discardedReads, stList *discardedAlignments);
bool downsampleBamChunkReadWithVcfEntrySubstringsViaFullReadLengthLikelihood(int64_t intendedDepth,
                                                                             stList *chunkVcfEntries,
                                                                             stList *inputReads,
                                                                             stList *maintainedReads,
                                                                             stList *discardedReads);

void writeHaplotaggedBam(BamChunk *bamChunk, char *inputBamLocation, char *outputBamFileBase,
                         stSet *readsInH1, stSet *readsInH2, Params *params, char *logIdentifier);

char *getSequenceFromReference(char *fastaFile, char *contig, int64_t startPos, int64_t endPosExcl);

int64_t getAlignedReadLength(bam1_t *aln);

int64_t getAlignedReadLength2(bam1_t *aln, int64_t *start_softclip, int64_t *end_softclip);

int64_t getAlignedReadLength3(bam1_t *aln, int64_t *start_softclip, int64_t *end_softclip, bool boundaryAtMatch);

void countIndels(uint32_t *cigar, uint32_t ncigar, int64_t *numInsertions, int64_t *numDeletions);


/*
 * Writes all supplemental information for a chunk
 */
void poa_writeSupplementalChunkInformation(char *outputBase, int64_t chunkIdx,
                                           BamChunk *bamChunk, Poa *poa, stList *reads, Params *params,
                                           bool outputPoaDOT, bool outputPoaCSV, bool outputRepeatCounts);
void poa_writeSupplementalChunkInformation2(char *outputBase, char *haplotypeIdentifier, int64_t chunkIdx,
                                            BamChunk *bamChunk, Poa *poa, stList *reads, Params *params,
                                            bool outputPoaDOT, bool outputPoaCSV, bool outputRepeatCounts);
void poa_writeSupplementalChunkInformationDiploid(char *outputBase, int64_t chunkIdx,
                                                  BamChunk *bamChunk, stGenomeFragment *genomeFragment, Poa *poaH1, Poa *poaH2, stList *bamChunkReads,
                                                  stSet *readsInHap1, stSet *readsInHap2, Params *params, bool outputPoaDOT, bool outputPoaCSV,
                                                  bool outputRepeatCounts, bool outputHaplotypedReadIdCsv, bool outputHaplotypedBam, char *logIdentifier);



#endif /* ST_RP_HMM_H_ */
