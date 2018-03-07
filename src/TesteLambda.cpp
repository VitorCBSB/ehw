//============================================================================
// Name        : TesteLambda.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <functional>
#include <vector>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdint.h>
#include <stdio.h>
#include <cassert>
#include "GenericGA.h"
#include "random.h"
#include "Utils.h"
#include "hps_0.h"

#define REG_BASE 0xFF200000
#define REG_SPAN 0x00200000

enum Function {
	AND,
	OR,
	NOT,
	XOR,
	XNOR,
	NAND,
	NOR
};

struct Cell {
	Function func;
	std::vector<unsigned int> inputs;
};

struct GeneticParams {
	unsigned int r, c, numIn, numOut, leNumIn;
};

struct Chromosome {
	std::vector<std::vector<Cell>> cells;
	std::vector<unsigned int> outputs;
};

struct CircuitAnalysis {
	int maxDepth;
	int logicGatesUsed;
	int transistorsUsed;
};

template <typename T, typename U>
std::tuple<T, U> tup(T t, U u) {
    return std::make_tuple(t, u);
}

bool dominates(CircuitAnalysis a, CircuitAnalysis b) {
	return (a.maxDepth < b.maxDepth || a.logicGatesUsed < b.logicGatesUsed || a.transistorsUsed < b.transistorsUsed)
			&& !(b.maxDepth < a.maxDepth || b.logicGatesUsed < a.logicGatesUsed || b.transistorsUsed < a.transistorsUsed);
}

CircuitAnalysis  mergeAnalysis(CircuitAnalysis a, CircuitAnalysis b) {
	CircuitAnalysis result;
	result.maxDepth = std::max(a.maxDepth, b.maxDepth);
    result.logicGatesUsed = a.logicGatesUsed + b.logicGatesUsed;
    result.transistorsUsed = a.transistorsUsed + b.transistorsUsed;
	return result;
}

Chromosome  makeChromosome
		( std::vector<std::vector<Cell>> cells
        , std::vector<unsigned int> outputs
		) {
	Chromosome res;
	res.cells = cells;
	res.outputs = outputs;
	return res;
}

// bool vector has at most 32 elements
uint32_t convertUnit(std::vector<bool> v) {
    uint32_t res = 0;
    for (unsigned int i = 0; i < v.size(); i++) {
        res |= (v[i] << i);
    }
    return res;
}

// bits range = [0..32]
std::vector<bool> convertToBits(uint32_t val, uint32_t bits) {
    std::vector<bool> res;
    for (unsigned int i = 0; i < bits; i++) {
        res.push_back((val & (1 << i)) >> i);
    }
    return res;
}

std::vector<uint32_t> convertToPacked(std::vector<bool> v) {
    return map(convertUnit, chunksOf(v, 32));
}

std::vector<bool> serializeCell(GeneticParams params, Cell cell) {
	std::vector<bool> result;

	auto numPinos = params.numIn + params.c * params.r;
	auto bitsPinos = (uint32_t) ceil(log2(numPinos));
	auto bitsFunc = 3;

	auto inputs = map([=](unsigned int in) {
		return convertToBits(in, bitsPinos);
	}, cell.inputs);

	auto func = convertToBits(cell.func, bitsFunc);

	for (auto in : inputs) {
		result.insert(result.end(), in.begin(), in.end());
	}
	result.insert(result.end(), func.begin(), func.end());

	return result;
}

std::vector<bool> rawSerialize(GeneticParams params, Chromosome chrom) {
 	std::vector<uint32_t> result;

	auto numPinos = params.numIn + params.c * params.r;
	auto bitsPinos = (uint32_t) ceil(log2(numPinos));

	std::vector<bool> totalBits;

	for (unsigned int j = 0; j < params.c; j++) {
		for (unsigned int i = 0; i < params.r; i++) {
			auto cellBits = serializeCell(params, chrom.cells[i][j]);
			totalBits.insert(totalBits.end(), cellBits.begin(), cellBits.end());
		}
	}

	auto outBits = map([=](unsigned int out) {
		return convertToBits(out, bitsPinos);
	}, chrom.outputs);

	for (auto out : outBits) {
        totalBits.insert(totalBits.end(), out.begin(), out.end());
	}

	return totalBits;
}

std::vector<uint32_t> serialize(GeneticParams params, Chromosome chrom) {
	return convertToPacked(rawSerialize(params, chrom));
}

Cell makeCell(Function func, std::vector<unsigned int> inputs) {
	Cell res;
	res.func = func;
	res.inputs = inputs;
	return res;
}

unsigned int transistors(Function func) {
	switch(func) {
	case AND: return 6;
	case OR: return 6;
	case XOR: return 6;
	case NOT: return 2;
	case NAND: return 4;
	case XNOR: return 8;
	case NOR: return 4;
	default: return 6; // Shutting the compiler up...
	}
}

uint32_t  firstBitOfEachOutput(std::vector<std::bitset<8>> answer) {
	uint32_t result = 0;
	for (unsigned int i = 0; i < answer.size(); i++) {
		result = result + (answer[i][0] << i);
	}
	return result;
}

CircuitAnalysis recursiveBuild
		( GeneticParams params
		, const Chromosome& chrom
        , std::vector<std::vector<bool>>& visitedNodes
		, std::vector<std::vector<CircuitAnalysis>>& analysisMatrix
		, unsigned int input
		) {
    CircuitAnalysis res;
    res.maxDepth = 0;
    res.logicGatesUsed = 0;
    res.transistorsUsed = 0;
	if (input < params.numIn) {
		return res;
	}

	unsigned int actualInput = input - params.numIn;
	unsigned int i = actualInput % params.r;
	unsigned int j = actualInput / params.r;

	if (visitedNodes[i][j]) {
		res.maxDepth = analysisMatrix[i][j].maxDepth;
		return res;
	}

	res.logicGatesUsed = 1;
	res.transistorsUsed = transistors(chrom.cells[i][j].func);

	if (chrom.cells[i][j].func == NOT) {
		res = mergeAnalysis
				( res
                , recursiveBuild
                    ( params
                    , chrom
					, visitedNodes
					, analysisMatrix
					, chrom.cells[i][j].inputs[0]
                    )
                );
	} else {
		res = fold(mergeAnalysis, res, map([&](unsigned int curInput) {
			return recursiveBuild(params, chrom, visitedNodes, analysisMatrix, curInput);
		}, chrom.cells[i][j].inputs));
	}
	res.maxDepth = res.maxDepth + 1;

	visitedNodes[i][j] = true;
	analysisMatrix[i][j] = res;

	return res;
}

CircuitAnalysis  circuitAnalysis(GeneticParams params, const Chromosome& chrom) {
	std::vector<std::vector<CircuitAnalysis>> analysisMatrix(params.r, std::vector<CircuitAnalysis>(params.c));
	std::vector<std::vector<bool>> visitedNodes(params.r, std::vector<bool>(params.c));

	CircuitAnalysis init;
	init.maxDepth = init.logicGatesUsed = init.transistorsUsed = 0;

	return fold(mergeAnalysis, init, map([&](unsigned int output) {
		return recursiveBuild(params, chrom, visitedNodes, analysisMatrix, output);
	}, chrom.outputs));
}

std::tuple<unsigned int, unsigned int> indexToCoordinate(unsigned int index, unsigned int r) {
	return std::make_tuple(index % r, index / r);
}

unsigned int rawCoordinateToIndex(std::tuple<unsigned int, unsigned int> coordinate, unsigned int r) {
	return std::get<0>(coordinate) + r * std::get<1>(coordinate);
}

unsigned int coordinateToIndex
    ( std::tuple<unsigned int, unsigned int> coordinate
    , unsigned int r
    , unsigned int numIn
    ) {
    return numIn + rawCoordinateToIndex(coordinate, r);
}

unsigned int fitInLargerIndex(unsigned int index, unsigned int oldR, unsigned int newR, unsigned int numIn, unsigned int lineOffset) {
	if (index >= numIn) {
		index -= numIn;
		auto coord = indexToCoordinate(index, oldR);
		std::get<0>(coord) += lineOffset;
		return coordinateToIndex(coord, newR, numIn);
	}
	return index;
}

Cell fitInLargerCell(Cell cell, unsigned int oldR, unsigned int newR, unsigned int numIn) {
	return makeCell(cell.func, map([=](unsigned int input) {
		return fitInLargerIndex(input, oldR, newR, numIn, 0);
	}, cell.inputs));
}

// ALL parameters in largerParams are assumed to be greater or equal to params.
Chromosome fitInLargerChrom(Chromosome chrom, unsigned int outputNum, GeneticParams params, GeneticParams largerParams) {
	Chromosome result;

	result.cells = replicate(largerParams.r, replicate(largerParams.c, makeCell(AND, replicate(params.leNumIn, (unsigned int) 0))));
	for (unsigned int i = 0; i < params.r; i++) {
		for (unsigned int j = 0; j < params.c; j++) {
			result.cells[i][j] = fitInLargerCell(chrom.cells[i][j], params.r, largerParams.r, largerParams.numIn);
		}
	}

	for (unsigned int i = 0; i < outputNum; i++) {
		result.outputs.push_back(0);
	}
	result.outputs.insert(result.outputs.end(), chrom.outputs.begin(), chrom.outputs.end());
	for (unsigned int i = outputNum; i < params.numOut + outputNum; i++) {
		result.outputs[i] = fitInLargerIndex(result.outputs[i], params.r, largerParams.r, largerParams.numIn, 0);
	}
	for (unsigned int i = 0; i < largerParams.numOut - (params.numOut + outputNum); i++) {
		result.outputs.push_back(0);
	}

	return result;
}

std::string showFunction(Function f) {
	switch(f) {
	case AND: return "AND";
	case OR: return "OR";
	case NOT: return "NOT";
	case XOR: return "XOR";
	case XNOR: return "XNOR";
	case NAND: return "NAND";
	case NOR: return "NOR";
	}
	return "What";
}

std::string showInt(unsigned int i) {
	std::ostringstream ss;
	ss << i;
	return ss.str();
}

std::string showInput(GeneticParams params, unsigned int i) {
	std::string s;
	if (i < params.numIn) {
		s += "#";
		s += showInt(i);
	} else {
		i -= params.numIn;
        auto col = i / params.r;
        auto row = i % params.r;
        s += "(";
        s += showInt(row);
        s += ", ";
        s += showInt(col);
        s += ")";
	}
	return s;
}

std::string showCell(GeneticParams params, Cell c) {
	auto s = showFunction(c.func);
	s += "[";
	for (unsigned int i = 0; i < c.inputs.size(); i++) {
		s += showInput(params, c.inputs[i]);
		if (i < c.inputs.size() - 1) {
			s += ", ";
		}
	}
	s += "]";
	return s;
}

std::vector<std::string> format(std::vector<std::string> vs) {
	auto sizes = map([=](std::string s) { return s.size(); }, vs);
	auto max = std::max_element(sizes.begin(), sizes.end());
	return map([=](std::string s) {
		auto origSize = s.size();
		for (unsigned int i = 0; i < *max - origSize; i++) {
			s += ' ';
		}
		return s;
	}, vs);
}

std::string showChromosome(GeneticParams params, Chromosome chrom) {
	std::string s;
	s += "Cells:\n";
	auto cellsS =
			map([=](std::vector<Cell> row) {
		return map([=](Cell c) { return showCell(params, c); }, row);
	}, chrom.cells);

	auto tCellsS = transpose(cellsS);

	auto tFormatted = map
			([=](std::vector<std::string> column) {
		return format(column);
	}, tCellsS);

	auto formatted = transpose(tFormatted);

	for (unsigned int i = 0; i < formatted.size(); i++) {
		for (unsigned int j = 0; j < formatted[0].size(); j++) {
			s += formatted[i][j];
			s += "   ";
		}
		s += '\n';
	}

	s += "Outputs:\n";
	for (unsigned int i = 0; i < chrom.outputs.size(); i++) {
		s += showInput(params, chrom.outputs[i]);
		s += ' ';
	}
	s += '\n';
	return s;
}

// This function takes the original GeneticParams
// used to reach each individual solution.
// The new one should be calculated after this function
// is called.
Chromosome  mergeChromosomes(GeneticParams params, std::vector<Chromosome> chroms) {
	auto zipped = zip(vectorFromTo(0, chroms.size()), chroms);
	auto newR = chroms.size() * params.r;

	auto transformed = map([=](std::tuple<unsigned int, Chromosome> tup) {
		auto transform = [=](unsigned int input) -> unsigned int {
			return fitInLargerIndex(input, params.r, newR, params.numIn, std::get<0>(tup) * params.r);
		};

		auto newCells = map([=](std::vector<Cell> cells) {
			return map([=](Cell cell) {
				return makeCell(cell.func, map(transform, cell.inputs));
			}, cells);
		}, std::get<1>(tup).cells);

		auto newOuts = map(transform, std::get<1>(tup).outputs);

		return std::make_tuple(std::get<0>(tup), makeChromosome(newCells, newOuts));
	}, zipped);

	std::vector<std::vector<Cell>> initCells(newR, std::vector<Cell>(params.c));
	std::vector<unsigned int> initOuts(params.numOut * chroms.size());

	return fold([=](Chromosome result, std::tuple<unsigned int, Chromosome> tup) {
		for (unsigned int i = 0; i < params.r; i++) {
			for (unsigned int j = 0; j < params.c; j++) {
				result.cells[i + std::get<0>(tup) * params.r][j] = std::get<1>(tup).cells[i][j];
			}
		}
		for (unsigned int i = 0; i < params.numOut; i++) {
            result.outputs[i + std::get<0>(tup) * params.numOut] = std::get<1>(tup).outputs[i];
		}
		return result;
	}, makeChromosome(initCells, initOuts), transformed);
}

std::vector<std::tuple<std::bitset<8>, std::bitset<8>, std::bitset<8>>>
    inputOutputValidSequences() {
    return map([](std::tuple<const char*, const char*, const char*> s)
            { return std::make_tuple(std::bitset<8>(std::get<0>(s)), std::bitset<8>(std::get<1>(s)), std::bitset<8>(std::get<2>(s))); },
            std::vector<std::tuple<const char*, const char*, const char*>>{
                std::make_tuple("00000010","00000000", "00000000"),
                std::make_tuple("00000001","00000001", "11111111"),
                std::make_tuple("00000010","00000000", "11111111"),
                std::make_tuple("00000011","00000000", "11111111"),
                std::make_tuple("00000000","00000001", "11111111"),
                std::make_tuple("00000011","00000001", "11111111"),
                std::make_tuple("00000010","00000001", "11111111"),
                std::make_tuple("00000011","00000000", "11111111"),
                std::make_tuple("00000010","00000000", "11111111"),
                std::make_tuple("00000000","00000001", "11111111"),
                std::make_tuple("00000010","00000000", "11111111"),
                std::make_tuple("00000001","00000001", "11111111"),
                std::make_tuple("00000000","00000001", "11111111"),
                std::make_tuple("00000001","00000000", "11111111"),
                std::make_tuple("00000000","00000000", "11111111"),
                std::make_tuple("00000001","00000000", "11111111")
            });
}

Function functionFromInt(unsigned int n) {
	switch (n) {
	case 0:
		return AND;
	case 1:
		return OR;
	case 2:
		return NOR;
	case 3:
		return XOR;
	case 4:
		return XNOR;
	case 5:
		return NAND;
	case 6:
		return NOR;
	default:
		return AND;
	}
}

std::function<bool(std::vector<bool>)>
		cellFunction(Function func) {
	switch (func) {
	case AND:
		return [](std::vector<bool> in) { return fold1(std::logical_and<bool>(), in); };
	case OR:
		return [](std::vector<bool> in) { return fold1(std::logical_or<bool>(), in); };
	case NOT:
		return [](std::vector<bool> in) { return !in[0]; };
	case XOR:
		return [](std::vector<bool> in) { return fold1(std::not_equal_to<bool>(), in); };
	case XNOR:
		return [](std::vector<bool> in) { return !fold1(std::not_equal_to<bool>(), in); };
	case NAND:
		return [](std::vector<bool> in) { return !fold1(std::logical_and<bool>(), in); };
	case NOR:
		return [](std::vector<bool> in) { return !fold1(std::logical_or<bool>(), in); };
	}
	// This is here to shut the compiler up. It should never actually be reached.
    return [](std::vector<bool> in) { return fold1(std::logical_and<bool>(), in); };
}

bool
	chooseInput
		( std::vector<std::vector<bool>> resultMatrix
		, std::bitset<8> circuitIn
        , unsigned int input
		, unsigned int numIn
		, unsigned int r
		) {
	if (input >= numIn) {
		auto actualInput = input - numIn;
		return resultMatrix[actualInput % r][actualInput / r];
	}
	return circuitIn[input];
}

std::vector<std::vector<bool>>
		initialResultMatrix(unsigned int r, unsigned int c) {
	std::vector<std::vector<bool>> result(r, std::vector<bool>(c));
	return result;
}

// I'm making the assumption that the matrix is
// feed-forward only.
std::vector<std::vector<bool>>
		evaluateMatrix(Chromosome chrom, std::bitset<8> in, unsigned int r, unsigned int c, unsigned int numIn) {
	auto resultMatrix = initialResultMatrix(r, c);

	for (unsigned int j = 0; j < c; j++) {
		for (unsigned int i = 0; i < r; i++) {
			auto inputs = map(
					[=](unsigned int in_) {
                        return chooseInput(resultMatrix, in, in_, numIn, r);
                    }, chrom.cells[i][j].inputs);
			resultMatrix[i][j] = cellFunction(chrom.cells[i][j].func)(inputs);
		}
	}

	return resultMatrix;
}

double
    evaluateChromosome
		( GeneticParams params
        , Chromosome chrom
		, std::vector<std::bitset<8>> expectedResult
		) {
    unsigned int distanceSum =
        fold([=](unsigned int totDistanceSum, unsigned int input) {
            auto resultMatrix = evaluateMatrix(chrom, std::bitset<8>(input), params.r, params.c, params.numIn);
            return totDistanceSum +
                fold([=](unsigned int distanceSum, unsigned int outIndex) {
                        return distanceSum
                                + (expectedResult[input][outIndex]
                                                 != chooseInput(resultMatrix, std::bitset<8>(input), chrom.outputs[outIndex], params.numIn, params.r));
                }, 0, vectorFromTo(0, chrom.outputs.size()));
        }, 0, vectorFromTo(0, pow(2, params.numIn)));

    if (distanceSum == 0) {
        return 2000000.0;
    }
    return 1.0 / (double) distanceSum;
}

std::function<double(Chromosome)>
		makeFitnessFunc
			( GeneticParams params
            , std::vector<std::bitset<8>> expectedResult
			) {
	return [=](Chromosome chrom) {
		return evaluateChromosome(params, chrom, expectedResult);
	};
}

// Side-effectful!! Writes to FPGA ports.
void sendVectorToFPGA(std::vector<uint32_t> vec, void* fpgaMemory) {
	std::vector<int> segmentAddrs = {
		CHROM_SEG_0_BASE, CHROM_SEG_1_BASE, CHROM_SEG_2_BASE, CHROM_SEG_3_BASE, CHROM_SEG_4_BASE,
		CHROM_SEG_5_BASE, CHROM_SEG_6_BASE, CHROM_SEG_7_BASE, CHROM_SEG_8_BASE, CHROM_SEG_9_BASE,
		CHROM_SEG_10_BASE, CHROM_SEG_11_BASE, CHROM_SEG_12_BASE, CHROM_SEG_13_BASE, CHROM_SEG_14_BASE,
		CHROM_SEG_15_BASE, CHROM_SEG_16_BASE, CHROM_SEG_17_BASE, CHROM_SEG_18_BASE, CHROM_SEG_19_BASE,
		CHROM_SEG_20_BASE, CHROM_SEG_21_BASE, CHROM_SEG_22_BASE, CHROM_SEG_23_BASE, CHROM_SEG_24_BASE,
		CHROM_SEG_25_BASE, CHROM_SEG_26_BASE, CHROM_SEG_27_BASE, CHROM_SEG_28_BASE, CHROM_SEG_29_BASE,
		CHROM_SEG_30_BASE
	};

	// Sending serialized chromosome to FPGA.
	for (unsigned int i = 0; i < vec.size(); i++) {
		void* segmentAddr = (uint8_t*) fpgaMemory + segmentAddrs[i];
		*(uint32_t*) segmentAddr = vec[i];
	}
}

void sendChromosomeToFPGA(Chromosome chromosome, GeneticParams params, void* fpgaMemory) {
	sendVectorToFPGA(serialize(params, chromosome), fpgaMemory);
}

uint32_t sendVectorAndGetErrorSum
    ( std::vector<uint32_t> vec
    , void* fpgaMemory
    ) {

    std::vector<int> errorSumAddrs = {
        ERROR_SUM_0_BASE, ERROR_SUM_1_BASE, ERROR_SUM_2_BASE, ERROR_SUM_3_BASE,
        ERROR_SUM_4_BASE, ERROR_SUM_5_BASE, ERROR_SUM_6_BASE, ERROR_SUM_7_BASE
    };

    void* doneProcessingFeedbackAddr = (uint8_t*) fpgaMemory + DONE_PROCESSING_FEEDBACK_BASE;
    *(uint32_t*) doneProcessingFeedbackAddr = 0;

    void* readyToProcessAddr = (uint8_t*) fpgaMemory + READY_TO_PROCESS_BASE;
    while ((*(uint32_t*) readyToProcessAddr) != 1);

    sendVectorToFPGA(vec, fpgaMemory);

    void* startProcessingAddr = (uint8_t*) fpgaMemory + START_PROCESSING_CHROM_BASE;
    *(uint32_t*) startProcessingAddr = 1;

    void* doneProcessingAddr = (uint8_t*) fpgaMemory + DONE_PROCESSING_CHROM_BASE;
    while ((*(uint32_t*) doneProcessingAddr) != 1) {
        // Esse print esta aqui porque o otimizador do g++ faz o programa
        // ficar preso num loop infinito se isso nao estiver aqui.
        std::cout << "";
    }

    uint32_t chromErrorSum = 0;
    for (auto addr : errorSumAddrs) {
        void* chromErrorSumAddr = (uint8_t*) fpgaMemory + addr;
        chromErrorSum += *(uint32_t*) chromErrorSumAddr;
    }

    // Esse print esta aqui porque o otimizador do g++ faz o programa
    // ficar preso num loop infinito se isso nao estiver aqui.
    std::cout << "";

    *(uint32_t*) startProcessingAddr = 0;
    *(uint32_t*) doneProcessingFeedbackAddr = 1;

    return chromErrorSum;
}

uint32_t sendChromosomeAndGetErrorSum
    ( Chromosome chrom
    , GeneticParams params
    , void* fpgaMemory
    ) {
    auto serial = serialize(params, chrom);
    return sendVectorAndGetErrorSum(serial, fpgaMemory);
}

std::function<double(Chromosome)>
		makeFPGAFitnessFunc
			( GeneticParams params
			, void* fpgaMemory
            ) {
	return [=](Chromosome chrom) {
	    auto sum = sendChromosomeAndGetErrorSum(chrom, params, fpgaMemory);

		if (sum == 0) {
            return 2000000.0;
		}
		return 1.0 / (double) sum;
	};
}

std::function<double(Chromosome)>
		makeGrowingFPGAFitnessFunc
			( GeneticParams params
            , GeneticParams largerParams
			, void* fpgaMemory
            ) {
	return [=](Chromosome chrom) {
	    auto largerChrom = fitInLargerChrom(chrom, 0, params, largerParams);
	    auto sum = sendChromosomeAndGetErrorSum(largerChrom, largerParams, fpgaMemory);

		if (sum == 0) {
            return 2000000.0;
		}
		return 1.0 / (double) sum;
	};
}

std::function<double(Chromosome)>
		makeFPGASeparateFitnessFunc
			( GeneticParams params
			, unsigned int outputNum
			, unsigned int numOutputs
			, void* fpgaMemory
            ) {
	return [=](Chromosome chrom) {
		auto largerParams = params;
		largerParams.r = params.r * numOutputs;
		largerParams.numOut = numOutputs;

		auto sum = sendChromosomeAndGetErrorSum(fitInLargerChrom(chrom, outputNum, params, largerParams), largerParams, fpgaMemory);

		if (sum == 0) {
            return 2000000.0;
		}
		return 1.0 / (double) sum;
	};
}

std::function<double(Chromosome)>
		makeOptimizingFitnessFunc
			( GeneticParams params
            , std::vector<std::bitset<8>> expectedResult
			) {
	return [=](Chromosome chrom) {
		auto chromScore = evaluateChromosome(params, chrom, expectedResult);
		if (chromScore != 2000000.0) {
			return 0.0;
		}

		auto analysis = circuitAnalysis(params, chrom);

		return (1.0 / analysis.maxDepth + 1.0 / analysis.logicGatesUsed + 1.0 / analysis.transistorsUsed);
	};
}

std::function<double(Chromosome)>
	makeOptimizingFPGAFitnessFunc
		( GeneticParams params
        , void* fpgaMemory
		) {
	return [=](Chromosome chrom) {
	    auto sum = sendChromosomeAndGetErrorSum(chrom, params, fpgaMemory);

		if (sum != 0) {
            return 0.0;
		}
		auto analysis = circuitAnalysis(params, chrom);

		return (1.0 / analysis.maxDepth + 1.0 / analysis.logicGatesUsed + 1.0 / analysis.transistorsUsed);
	};
}

RNGFUNC(Function) randomFunc() {
	return rmap<random_type, Function>([](random_type r) {
		return functionFromInt(r % 7);
	}, getRandom());
}

RNGFUNC(unsigned int)  randomOutput(GeneticParams params) {
	return bind
			( getRandom()
            , [=](random_type rand) {
		return pure(rand % (params.numIn + params.r * params.c));
	});
}

RNGFUNC(std::vector<unsigned int>)  mutateOutput
		( std::vector<unsigned int> outputs
        , unsigned int pointToMutate
        , GeneticParams params) {
	return bind
			( randomOutput(params)
            , [=](unsigned int newOut) mutable {
		outputs[pointToMutate] = newOut;
		return pure(outputs);
	});
}

RNGFUNC(Cell)  randomCell
		( GeneticParams params
		) {
	return bind(randomFunc(), [=](Function randFunc) {
		return bind(sequence(replicate(params.leNumIn, randomOutput(params)))
				, [=](std::vector<unsigned int> randomInputs) {
			return pure(makeCell(randFunc, randomInputs));
		});
	});
}

RNGFUNC(std::vector<std::vector<Cell>>)  mutateGrid
		( std::vector<std::vector<Cell>> grid
        , random_type pointToMutate
		, GeneticParams params
		) {
	return bind(getRandom(), [=](random_type rand) mutable {
            auto attrToMutate = rand % (params.leNumIn + 1);
            auto col = pointToMutate / params.r;
            auto row = pointToMutate % params.r;
            if (attrToMutate == 0) {
            	return bind(randomFunc(), [=](Function rFunc) mutable {
            		grid[row][col].func = rFunc;
            		return pure(grid);
            	});
            } else {
                return bind(randomOutput(params), [=](unsigned int rIn) mutable {
                    grid[row][col].inputs[(attrToMutate - 1)] = rIn;
                    return pure(grid);
                });
            }
        });
}

std::function<RNGFUNC(Chromosome)(Chromosome)>
	makeMutation(GeneticParams params, float mutationPercent) {
    auto totalElements = params.r * params.c * (params.leNumIn + 1) + params.numOut;
    auto elementsToMutate = std::ceil(totalElements * mutationPercent);

	return [=](Chromosome chrom) {
		return bind(sequence(replicate(elementsToMutate, getRandom())),
				[=](std::vector<random_type> rands) mutable {
			return foldM([=](Chromosome c, random_type rand) mutable {
                auto pointToMutate = rand % (params.c * params.r + params.numOut);
                if (pointToMutate < params.c * params.r) {
                    return bind
                            ( mutateGrid(chrom.cells, pointToMutate, params)
                            , [=](std::vector<std::vector<Cell>> newGrid) mutable {
                        c.cells = newGrid;
                        return pure(c);
                    });
                } else {
                    return bind
                            ( mutateOutput(chrom.outputs, pointToMutate - (params.c * params.r), params)
                            , [=](std::vector<unsigned int> newOuts) mutable {
                        c.outputs = newOuts;
                        return pure(c);
                    });
                }
			}, chrom, rands);
		});
	};
}

RNGFUNC(std::vector<Cell>)  randomColumn(GeneticParams params) {
	return sequence(replicate(params.r, randomCell(params)));
}

RNGFUNC(std::vector<std::vector<Cell>>)  randomCells(GeneticParams params) {
	return bind(mapM([=](unsigned int unused) {
		return randomColumn(params);
	}, vectorFromTo(0, params.c)), [=](std::vector<std::vector<Cell>> grid) {
		return pure(transpose(grid));
	});
}

RNGFUNC(std::vector<unsigned int>)  randomOutputs(GeneticParams params) {
	return sequence(replicate(params.numOut, randomOutput(params)));
}

RNGFUNC(Chromosome)  randomChrom(GeneticParams params) {
	return bind(randomOutputs(params), [=](std::vector<unsigned int> randomOuts) {
		return bind(randomCells(params), [=](std::vector<std::vector<Cell>> randomCs) {
			return pure(makeChromosome(randomCs, randomOuts));
		});
	});
}

template <typename F>
std::function<bool(GAState<Evaluated<Chromosome>>)>
    makeCorrectTermination(F fitness) {

    static_assert(std::is_convertible<F, std::function<double(Chromosome)>> ::value,
                "any's function must be of type Chromosome -> double");

    return [=](GAState<Evaluated<Chromosome>> state) {
        if (state.population[0].score >= 2000000) {
            auto verifiedScore = fitness(state.population[0].value);
            if (verifiedScore >= 2000000) {
                printf("%d %g\n", state.generation, state.population[0].score);
                return false;
            }
            std::cout << "Failed re-verification." << std::endl;
            return true;
        }
        if (state.generation % 100 == 0) {
            printf("%d %g\n", state.generation, state.population[0].score);
        }
        return state.generation < 200000;
    };
}

bool correctTermination(GAState<Evaluated<Chromosome>> state) {
    if (state.generation % 10 == 0 || state.population[0].score >= 2000000) {
        printf("%d %g\n", state.generation, state.population[0].score);
    }
	return state.generation < 50000 && state.population[0].score < 2000000;
}

bool  optimizingTermination(GAState<Evaluated<Chromosome>> state) {
	printf("%d %g\n", state.generation, state.population[0].score);
	return state.generation < 50000;
}

template <typename F>
RNGFUNC(GAState<Evaluated<Chromosome>>)
	GARoutineWithInitial
		( GeneticParams params
		, F fitnessFunc
		, Chromosome initial
        ) {

    static_assert(std::is_convertible<F, std::function<double(Chromosome)>> ::value,
                "fpgaGARoutine's fitness function must be of type T -> bool");

    auto mutationFunc = makeMutation(params, 0.15);
    auto strategy = lambdaPlusN<Chromosome>(fitnessFunc, mutationFunc, 4);
    auto gaFunc = makeGAFunction<Evaluated<Chromosome>>(strategy);
    auto termination = makeCorrectTermination(fitnessFunc);

    GAState<Evaluated<Chromosome>> init;
    init.generation = 0;
    init.population = { makeEvaluated(initial, fitnessFunc(initial)) };

    return iterateWhileM(termination, gaFunc, init);
}

RNGFUNC(std::vector<GAState<Evaluated<Chromosome>>>)
	simulatedGARoutine
		( GeneticParams params
        , std::vector<std::vector<std::bitset<8>>> answers
		) {
    return mapM([=](std::vector<std::bitset<8>> answer) {
        return bind
        		( randomChrom(params)
                , [=](Chromosome randomChromosome) {
            return GARoutineWithInitial
                       (params
                       , makeFitnessFunc(params, answer)
                       , randomChromosome
                       );
        });
	}, answers);
}

RNGFUNC(std::vector<GAState<Evaluated<Chromosome>>>)
	fpgaSeparateGARoutine
		( GeneticParams params
		, unsigned int totalNumOutputs
        , void* fpgaMemory) {
	return mapM([=](unsigned int outputNum) {
        return bind
        		( randomChrom(params)
                , [=](Chromosome randomChromosome) {
            return GARoutineWithInitial
                       (params
                       , makeFPGASeparateFitnessFunc(params, outputNum, totalNumOutputs, fpgaMemory)
                       , randomChromosome
                       );
        });
	}, vectorFromTo(0, totalNumOutputs));
}

RNGFUNC(GAState<Evaluated<Chromosome>>)
	fpgaGARoutine
		( GeneticParams params
        , void* fpgaMemory
        ) {
    return bind
            ( randomChrom(params)
            , [=](Chromosome randomChromosome) {
        return GARoutineWithInitial
                   ( params
                   , makeFPGAFitnessFunc(params, fpgaMemory)
                   , randomChromosome
                   );
    });
}

RNGFUNC(GAState<Evaluated<Chromosome>>)
	fpgaGrowingGARoutine
		( GeneticParams params
        , GeneticParams largerParams
        , void* fpgaMemory
        ) {
    return bind
            ( randomChrom(params)
            , [=](Chromosome randomChromosome) {
        return GARoutineWithInitial
                   ( params
                   , makeGrowingFPGAFitnessFunc(params, largerParams, fpgaMemory)
                   , randomChromosome
                   );
    });
}

void* openFPGAMemory() {
	int fd = open("/dev/mem", O_RDWR | O_SYNC);
	if (!fd) {
		printf("Could not open /dev/mem.\n");
		exit(1);
	}

	void* virtualBase = mmap(NULL, REG_SPAN, PROT_READ | PROT_WRITE, MAP_SHARED, fd, REG_BASE);
	if (!virtualBase) {
		printf("Could not map physical memory into virtual.");
		exit(1);
	}
	return virtualBase;
}

Chromosome fourWayAndSolution(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(AND, { 0, 1 });
	result.cells[1][0] = makeCell(AND, { 2, 3 });
	result.cells[0][1] = makeCell(AND,
			{ rawCoordinateToIndex(std::make_tuple(0, 0), r) + numIn
					, rawCoordinateToIndex(std::make_tuple(1, 0), r) + numIn });
	result.cells[0][2] = makeCell(NOT,
			{ rawCoordinateToIndex(std::make_tuple(0, 1), r) + numIn, 0 });

	result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(0, 2), r) + numIn);

	return result;
}

Chromosome twoAndSolution(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(AND, { 0, 1 });
	result.cells[0][1] = makeCell(NOT,
			{ rawCoordinateToIndex(std::make_tuple(0, 0), r) + numIn
					, 0 });

	result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(0, 0), r) + numIn);

	return result;
}

Chromosome SRLatch(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ rawCoordinateToIndex(std::make_tuple(1, 0), r) + numIn
					, 0 });
	result.cells[1][0] = makeCell(NOR,
			{ rawCoordinateToIndex(std::make_tuple(0, 0), r) + numIn
					, 1 });

	result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(1, 0), r) + numIn);
	result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(0, 0), r) + numIn);

	return result;
}

Chromosome oscillator(unsigned int r, unsigned int c, unsigned int numIn) {
 	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ rawCoordinateToIndex(std::make_tuple(0, 0), r) + numIn
					, coordinateToIndex(tup(0, 0), r, numIn) });
	result.cells[1][0] = makeCell(NAND,
			{ rawCoordinateToIndex(std::make_tuple(1, 0), r) + numIn
					, coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[0][1] = makeCell(NOT,
			{ rawCoordinateToIndex(std::make_tuple(0, 1), r) + numIn
					, 1 });
	result.cells[1][1] = makeCell(NOR,
			{ rawCoordinateToIndex(std::make_tuple(1, 1), r) + numIn
					, coordinateToIndex(tup(1, 1), r, numIn) });

	result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(0, 0), r) + numIn);
	result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(1, 0), r) + numIn);
	result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(0, 1), r) + numIn);
	result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(1, 1), r) + numIn);

	return result;
}

Chromosome Circuito1(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ 1, coordinateToIndex(std::make_tuple(1, 1), r, numIn) });
	result.cells[1][0] = makeCell(XNOR,
			{ rawCoordinateToIndex(std::make_tuple(0, 0), r) + numIn
					, 1 });
	result.cells[0][1] = makeCell(NOR, { 0, coordinateToIndex(std::make_tuple(0, 0), r, numIn) });
	result.cells[1][1] = makeCell(AND, { coordinateToIndex(std::make_tuple(0, 1), r, numIn),
	                                     coordinateToIndex(std::make_tuple(1, 0), r, numIn)});

	result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(0, 0), r) + numIn);
	result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(0, 1), r) + numIn);

	return result;
}

Chromosome Circuito2(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][0] = makeCell(NOR,
			{ 1, coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[0][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][1] = makeCell(NOR,
	        { 0, coordinateToIndex(std::make_tuple(1, 0), r, numIn)});

	result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(0, 0), r) + numIn);
	result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(1, 1), r) + numIn);
	//result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(0, 1), r) + numIn);
	//result.outputs.push_back(rawCoordinateToIndex(std::make_tuple(1, 0), r) + numIn);

	return result;
}

Chromosome Circuito3(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ coordinateToIndex(tup(1, 0), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][0] = makeCell(XOR,
			{ 0, coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[0][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), 1 });
	result.cells[1][1] = makeCell(AND,
	        { coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(0, 1), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(0, 1), r, numIn));
	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));

	return result;
}

Chromosome Circuito4(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ 1, coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][0] = makeCell(OR,
			{ 0, coordinateToIndex(tup(0, 0), r, numIn) });
	result.cells[0][1] = makeCell(XOR,
	        { coordinateToIndex(tup(0, 1), r, numIn), 1 });
	result.cells[1][1] = makeCell(NOR,
	        { coordinateToIndex(tup(1, 0), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(0, 1), r, numIn));
	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));

	return result;
}

Chromosome Circuito5(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(AND,
			{ coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][0] = makeCell(NOR,
			{ coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[0][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 1), r, numIn), 0 });
	result.cells[1][1] = makeCell(NOR,
	        { 1, coordinateToIndex(tup(0, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 1), r, numIn));
	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));

	return result;
}

Chromosome Circuito6(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ coordinateToIndex(tup(1, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][0] = makeCell(NOR,
			{ coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[0][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 1), r, numIn), 0 });
	result.cells[1][1] = makeCell(NOR,
	        { 1, coordinateToIndex(tup(0, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 1), r, numIn));
	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));

	return result;
}

Chromosome Circuito9(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ 1, coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][0] = makeCell(NOR,
			{ coordinateToIndex(tup(0, 0), r, numIn), 0 });
	result.cells[0][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(0, 0), r, numIn) });
	result.cells[1][1] = makeCell(NOR,
	        { 1, coordinateToIndex(tup(0, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));
	result.outputs.push_back(coordinateToIndex(tup(0, 1), r, numIn));

	return result;
}

Chromosome Circuito10(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(AND,
			{ coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][0] = makeCell(OR,
			{ coordinateToIndex(tup(0, 1), r, numIn), 1 });
	result.cells[0][1] = makeCell(NOR,
	        { 0, coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][1] = makeCell(XNOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 1), r, numIn));
	result.outputs.push_back(coordinateToIndex(tup(0, 1), r, numIn));

	return result;
}

Chromosome Circuito12(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(XOR,
			{ coordinateToIndex(tup(1, 0), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][0] = makeCell(OR,
			{ coordinateToIndex(tup(0, 1), r, numIn), 0 });
	result.cells[0][1] = makeCell(XNOR,
	        { coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), 1 });

	result.outputs.push_back(coordinateToIndex(tup(1, 1), r, numIn));
	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));

	return result;
}

Chromosome Circuito13(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ coordinateToIndex(tup(1, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][0] = makeCell(AND,
			{ coordinateToIndex(tup(0, 1), r, numIn), 1 });
	result.cells[0][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 1), r, numIn), 0 });
	result.cells[1][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), 0 });

	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));
	result.outputs.push_back(coordinateToIndex(tup(0, 1), r, numIn));

	return result;
}

Chromosome Circuito17(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ coordinateToIndex(tup(0, 0), r, numIn), coordinateToIndex(tup(0, 1), r, numIn) });
	result.cells[1][0] = makeCell(OR,
			{ coordinateToIndex(tup(1, 1), r, numIn), 0 });
	result.cells[0][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 1), r, numIn), 1 });

	result.outputs.push_back(coordinateToIndex(tup(1, 1), r, numIn));
	result.outputs.push_back(coordinateToIndex(tup(0, 1), r, numIn));

	return result;
}

Chromosome Circuito22(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(XOR,
			{ coordinateToIndex(tup(1, 1), r, numIn), 1 });
	result.cells[1][0] = makeCell(NOR,
			{ coordinateToIndex(tup(1, 1), r, numIn), coordinateToIndex(tup(0, 1), r, numIn) });
	result.cells[0][1] = makeCell(XOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), 1 });
	result.cells[1][1] = makeCell(NOR,
	        { coordinateToIndex(tup(1, 0), r, numIn), 0 });

	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));

	return result;
}

Chromosome Circuito23(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ coordinateToIndex(tup(1, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][0] = makeCell(XOR,
			{ coordinateToIndex(tup(0, 1), r, numIn), 1 });
	result.cells[0][1] = makeCell(XOR,
	        { coordinateToIndex(tup(1, 1), r, numIn), 1 });
	result.cells[1][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), 0 });

	result.outputs.push_back(coordinateToIndex(tup(1, 0), r, numIn));

	return result;
}

Chromosome Circuito24(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(XOR,
			{ coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][0] = makeCell(NAND,
			{ coordinateToIndex(tup(1, 1), r, numIn), 0 });
	result.cells[0][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), 0 });
	result.cells[1][1] = makeCell(NAND,
	        { 0, 1 });

	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));

	return result;
}

Chromosome Circuito25(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NAND,
			{ coordinateToIndex(tup(1, 0), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][0] = makeCell(XNOR,
			{ coordinateToIndex(tup(1, 1), r, numIn), coordinateToIndex(tup(0, 1), r, numIn) });
	result.cells[0][1] = makeCell(NAND,
	        { 0, 1 });
	result.cells[1][1] = makeCell(NOR,
	        { 0, coordinateToIndex(tup(0, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 0), r, numIn));

	return result;
}

Chromosome Circuito26(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ 0, coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][0] = makeCell(AND,
			{ coordinateToIndex(tup(0, 1), r, numIn), 0 });
	result.cells[0][1] = makeCell(XOR,
	        { 0, 1 });
	result.cells[1][1] = makeCell(NOR,
	        { coordinateToIndex(tup(1, 0), r, numIn) , coordinateToIndex(tup(0, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 1), r, numIn));

	return result;
}

Chromosome Circuito27(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(XNOR,
			{ coordinateToIndex(tup(1, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][0] = makeCell(OR,
			{ coordinateToIndex(tup(0, 1), r, numIn), 0 });
	result.cells[0][1] = makeCell(NOR,
	        { 0, coordinateToIndex(tup(0, 0), r, numIn) });
	result.cells[1][1] = makeCell(AND,
	        { 0, 1 });

	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));

	return result;
}

Chromosome Circuito28(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ coordinateToIndex(tup(1, 1), r, numIn), 1 });
	result.cells[1][0] = makeCell(NOR,
			{ coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[0][1] = makeCell(AND,
	        { 0, coordinateToIndex(tup(0, 0), r, numIn) });
	result.cells[1][1] = makeCell(NOR,
	        { 0, coordinateToIndex(tup(1, 1), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 0), r, numIn));

	return result;
}

Chromosome Circuito29(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(OR,
			{ coordinateToIndex(tup(1, 1), r, numIn), 1 });
	result.cells[1][0] = makeCell(XOR,
			{ coordinateToIndex(tup(0, 0), r, numIn), coordinateToIndex(tup(0, 1), r, numIn) });
	result.cells[0][1] = makeCell(NOR,
	        { 0, coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][1] = makeCell(XNOR,
	        { 0, 1 });

	result.outputs.push_back(coordinateToIndex(tup(1, 0), r, numIn));

	return result;
}

Chromosome Circuito30(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(XOR,
			{ coordinateToIndex(tup(1, 0), r, numIn), 1 });
	result.cells[1][0] = makeCell(NOR,
			{ 0, coordinateToIndex(tup(0, 1), r, numIn) });
	result.cells[0][1] = makeCell(NOR,
	        { coordinateToIndex(tup(1, 1), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][1] = makeCell(XOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), 1 });

	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));

	return result;
}

Chromosome Circuito31(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(XOR,
			{ coordinateToIndex(tup(1, 0), r, numIn), 1 });
	result.cells[1][0] = makeCell(NOR,
			{ 0, coordinateToIndex(tup(0, 1), r, numIn) });
	result.cells[0][1] = makeCell(NOR,
	        { coordinateToIndex(tup(1, 0), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][1] = makeCell(XOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), 1 });

	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));

	return result;
}

Chromosome Circuito32(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NOR,
			{ 0, 0 });
	result.cells[1][0] = makeCell(XNOR,
			{ coordinateToIndex(tup(1, 1), r, numIn), coordinateToIndex(tup(0, 1), r, numIn) });
	result.cells[0][1] = makeCell(NAND,
	        { 0, 1 });
	result.cells[1][1] = makeCell(AND,
	        { coordinateToIndex(tup(0, 0), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 0), r, numIn));

	return result;
}

Chromosome Circuito34(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(XNOR,
			{ coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][0] = makeCell(NAND,
			{ 0, coordinateToIndex(tup(0, 1), r, numIn) });
	result.cells[0][1] = makeCell(NAND,
	        { coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][1] = makeCell(XNOR,
	        { 0, 1 });

	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));

	return result;
}

Chromosome Circuito35(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(XOR,
			{ coordinateToIndex(tup(1, 0), r, numIn), 1 });
	result.cells[1][0] = makeCell(XOR,
			{ 1, coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[0][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), coordinateToIndex(tup(0, 0), r, numIn) });
	result.cells[1][1] = makeCell(NOR,
	        { 0, coordinateToIndex(tup(0, 1), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 0), r, numIn));

	return result;
}

Chromosome Circuito40(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NAND,
			{ coordinateToIndex(tup(0, 1), r, numIn), 1 });
	result.cells[1][0] = makeCell(NAND,
			{ coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(0, 0), r, numIn) });
	result.cells[0][1] = makeCell(XNOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][1] = makeCell(OR,
	        { 0, coordinateToIndex(tup(0, 1), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 0), r, numIn));

	return result;
}

Chromosome Circuito41(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NAND,
			{ coordinateToIndex(tup(1, 0), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][0] = makeCell(XNOR,
			{ coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[0][1] = makeCell(NOR,
	        { 0, coordinateToIndex(tup(0, 0), r, numIn) });
	result.cells[1][1] = makeCell(NAND,
	        { 1, coordinateToIndex(tup(0, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 0), r, numIn));

	return result;
}

Chromosome Circuito42(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(OR,
			{ coordinateToIndex(tup(1, 1), r, numIn), 1 });
	result.cells[1][0] = makeCell(NAND,
			{ coordinateToIndex(tup(0, 0), r, numIn), coordinateToIndex(tup(0, 1), r, numIn) });
	result.cells[0][1] = makeCell(NAND,
	        { 0, 1 });
	result.cells[1][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 0), r, numIn));

	return result;
}

Chromosome Circuito43(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(XNOR,
			{ coordinateToIndex(tup(1, 1), r, numIn), coordinateToIndex(tup(0, 1), r, numIn) });
	result.cells[1][0] = makeCell(NOR,
			{ coordinateToIndex(tup(0, 0), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[0][1] = makeCell(AND,
	        { 0, 1 });
	result.cells[1][1] = makeCell(NOR,
	        { coordinateToIndex(tup(1, 0), r, numIn), 1 });

	result.outputs.push_back(coordinateToIndex(tup(0, 0), r, numIn));

	return result;
}

Chromosome Circuito44(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NAND,
			{ 0, 1 });
	result.cells[1][0] = makeCell(OR,
			{ coordinateToIndex(tup(0, 1), r, numIn), 1 });
	result.cells[0][1] = makeCell(AND,
	        { coordinateToIndex(tup(0, 0), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });
	result.cells[1][1] = makeCell(NAND,
	        { coordinateToIndex(tup(0, 0), r, numIn), coordinateToIndex(tup(0, 1), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 1), r, numIn));

	return result;
}

Chromosome Circuito45(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(XNOR,
			{ coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][0] = makeCell(NOR,
			{ coordinateToIndex(tup(0, 0), r, numIn), 1 });
	result.cells[0][1] = makeCell(AND,
	        { 0, 1 });
	result.cells[1][1] = makeCell(OR,
	        { coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(1, 1), r, numIn));

	return result;
}

Chromosome DLatch(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(NAND,
			{ 1, 0 });
	result.cells[1][0] = makeCell(NAND,
			{ coordinateToIndex(tup(0, 0), r, numIn), 0 });
	result.cells[0][1] = makeCell(NAND,
	        { coordinateToIndex(tup(1, 1), r, numIn), coordinateToIndex(tup(0, 0), r, numIn) });
	result.cells[1][1] = makeCell(NAND,
	        { coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(0, 1), r, numIn));

	return result;
}

Chromosome JKLatch(unsigned int r, unsigned int c, unsigned int numIn) {
	Chromosome result;

	result.cells = replicate(r, replicate(c, makeCell(AND, replicate(2, (unsigned int) 0))));
	result.cells[0][0] = makeCell(AND,
			{ 0, coordinateToIndex(tup(0, 1), r, numIn) });
	result.cells[1][0] = makeCell(AND,
			{ coordinateToIndex(tup(1, 1), r, numIn), 1 });
	result.cells[0][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 0), r, numIn), coordinateToIndex(tup(1, 1), r, numIn) });
	result.cells[1][1] = makeCell(NOR,
	        { coordinateToIndex(tup(0, 1), r, numIn), coordinateToIndex(tup(1, 0), r, numIn) });

	result.outputs.push_back(coordinateToIndex(tup(0, 1), r, numIn));
	result.outputs.push_back(coordinateToIndex(tup(1, 1), r, numIn));

	return result;
}

std::vector<GeneticParams> growingGenParamVec(GeneticParams baseParams, unsigned int initialR, unsigned int finalR) {
    return map([=](unsigned int cur) {
        GeneticParams p = baseParams;
        p.r = cur;
        return p;
    }, vectorFromTo(initialR, finalR + 1));
}

int main() {
	GeneticParams params;
	params.r = 1; // Para cada solução individual.
	params.c = 1;
	params.numIn = 2 + 0; // + 1 para o clock.
	params.numOut = 1; // Para cada solução individual.
	params.leNumIn = 2;

	srand(time(NULL));
	auto initialRng = rand();

	printf("Seed: %d\n", initialRng);

	auto fpgaMem = openFPGAMemory();

	// Send input and output sequences to the FPGA.
	std::vector<int> outputSequencePorts = { EXPECTED_OUTPUT_0_BASE, EXPECTED_OUTPUT_1_BASE, EXPECTED_OUTPUT_2_BASE, EXPECTED_OUTPUT_3_BASE
				 };
	std::vector<int> inputSequencePorts = { INPUT_SEQUENCE_0_BASE, INPUT_SEQUENCE_1_BASE, INPUT_SEQUENCE_2_BASE, INPUT_SEQUENCE_3_BASE };
	std::vector<int> validOutputPorts = { VALID_OUTPUT_0_BASE, VALID_OUTPUT_1_BASE, VALID_OUTPUT_2_BASE, VALID_OUTPUT_3_BASE };

	auto io = inputOutputValidSequences();
	for (unsigned int i = 0; i < outputSequencePorts.size(); i++) {
	    void* inputAddress = (uint8_t*) fpgaMem + inputSequencePorts[i];
	    void* outputAddress = (uint8_t*) fpgaMem + outputSequencePorts[i];
	    void* validOutputAddress = (uint8_t*) fpgaMem + validOutputPorts[i];
	    auto curIndex = i * 4;
	    uint32_t inputVal =
	            std::get<0>(io[curIndex]).to_ulong()
	            | (std::get<0>(io[curIndex + 1]).to_ulong() << 8)
	            | (std::get<0>(io[curIndex + 2]).to_ulong() << 16)
	            | (std::get<0>(io[curIndex + 3]).to_ulong() << 24);
	    uint32_t outputVal =
	            std::get<1>(io[curIndex]).to_ulong()
	            | (std::get<1>(io[curIndex + 1]).to_ulong() << 8)
	            | (std::get<1>(io[curIndex + 2]).to_ulong() << 16)
	            | (std::get<1>(io[curIndex + 3]).to_ulong() << 24);
	    uint32_t validOutputVal =
	            std::get<2>(io[curIndex]).to_ulong()
	            | (std::get<2>(io[curIndex + 1]).to_ulong() << 8)
	            | (std::get<2>(io[curIndex + 2]).to_ulong() << 16)
	            | (std::get<2>(io[curIndex + 3]).to_ulong() << 24);
	    *(uint32_t*) inputAddress = inputVal;
	    *(uint32_t*) outputAddress = outputVal;
	    *(uint32_t*) validOutputAddress = validOutputVal;
	}

	/* Send a raw chromosome and evaluate once

	std::string bitstring = "001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000101001010000011000000000111010000100111110000110001011000001001100000011100000";
	std::reverse(bitstring.begin(), bitstring.end());
	auto data = std::vector<char>(bitstring.begin(), bitstring.end());
	auto fun = [](char c) {
	    return c == '1' ? true : false;
	};
	std::cout << sendVectorAndGetErrorSum(convertToPacked(map(fun, data)), fpgaMem) << std::endl;

	return 0;
	*/

	GeneticParams finalParams = params;
	finalParams.r = 25;

	auto b2c = [](bool b) {
	        return b ? '1' : '0';
	    };

	auto solution = whileFoldM([=](GeneticParams currentParams, GAState<Evaluated<Chromosome>> finalSolution) {
	    if (finalSolution.population[0].score < 2000000) {
	        std::cout << "Increasing number of rows by one, new value: " << currentParams.r + 1 << std::endl;
	        return true;
	    }

        printf("Solution:\n");
        printf("%s\n", showChromosome(currentParams, finalSolution.population[0].value).c_str());
        auto s = map(b2c, rawSerialize(currentParams, finalSolution.population[0].value));
        std::reverse(s.begin(), s.end());
        std::cout << "Normal bitstring:" << std::endl;
        std::cout << std::string(s.begin(), s.end()) << std::endl << std::endl;

        auto largerChrom = fitInLargerChrom(finalSolution.population[0].value, 0, currentParams, finalParams);
        s = map(b2c, rawSerialize(finalParams, largerChrom));
        std::reverse(s.begin(), s.end());
        std::cout << "Larger bitstring:" << std::endl;
        std::cout << std::string(s.begin(), s.end()) << std::endl << std::endl;

        auto fit = makeGrowingFPGAFitnessFunc(currentParams, finalParams, fpgaMem)(finalSolution.population[0].value);
        printf("Fitness recalculated: %g\n", fit);

	    return false;
	},
	[fpgaMem, finalParams](GeneticParams params) {
	    return fpgaGrowingGARoutine(params, finalParams, fpgaMem);
	}, growingGenParamVec(params, 3, finalParams.r));

	//auto solutions = fpgaGARoutineWithInitial(params, Circuito45(params.r, params.c, params.numIn), fpgaMem);
	//auto solutions = fpgaGARoutine(params, fpgaMem);

	evalState(solution, initialRng);

	/*
	unsigned int totalNumOut = 1;
	auto newParams = params;
	newParams.r *= totalNumOut;
	newParams.numOut = totalNumOut;

	auto optimizedSolution =
		bind
		( solutions
        , [=](std::vector<GAState<Evaluated<Chromosome>>> correctSolutions) {
            auto rawChroms = map([](GAState<Evaluated<Chromosome>> finalState) {
                return finalState.population[0].value;
            }, correctSolutions);

            auto merged = mergeChromosomes(params, rawChroms);

            auto mutationFunc = makeMutation(newParams, 0.05);
            auto fitnessFunc = makeOptimizingFPGAFitnessFunc(newParams, fpgaMem);
            auto strategy = lambdaPlusN<Chromosome>(fitnessFunc, mutationFunc, 4);
            auto gaFunc = makeGAFunction<Evaluated<Chromosome>>(strategy);

            GAState<Evaluated<Chromosome>> init;
            init.generation = 0;
            init.population = { makeEvaluated(merged, fitnessFunc(merged)) };

            return iterateWhileM(optimizingTermination, gaFunc, init);
	});
	*/

    /*
	auto analysis = circuitAnalysis(newParams, final.population[0].value);
	printf("Circuit analysis:\n");
	printf("Max depth: %d\n", analysis.maxDepth);
	printf("Logic gates: %d\n", analysis.logicGatesUsed);
	printf("Transistors: %d\n", analysis.transistorsUsed);
	*/

	return 0;
}
