#include "OT.h"
#include "OEP.h"
#include "PRNG.h"

// "On Arbitrary Waksman Networks and their Vulnerability"

// return sum(ceil(log2(i))) for i=1 to N
short gateSizeMap[] = { 0,1,3,5,8,11,14,17,21,25,29,33,37,41,45,49,54,59,64,
69,74,79,84,89,94,99,104,109,114,119,124,129,135,141,147,153,159,
165,171,177,183,189,195,201,207,213,219,225,231,237,243,249,255,
261,267,273,279,285,291,297,303,309,315,321,328,335,342,349,356,
363,370,377,384,391,398,405,412,419,426,433,440,447,454,461,468,
475,482,489,496,503,510,517,524,531,538,545,552,559,566 };

int getGateSize(int N)
{
    if (N < 100)
        return gateSizeMap[N - 1];
    int power = log2_32((uint32_t)N) + 1;
    return power * N + 1 - (1 << power);
}

void permutationToBits(int* permuIndices, int size, bool* bits)
{
    if (size == 2)
        bits[0] = permuIndices[0];
    if (size <= 2)
        return;

    int* invPermuIndices = new int[size];
    for (int i = 0; i < size; i++)
        invPermuIndices[permuIndices[i]] = i;

    bool odd = size & 1;

    // Solve the edge coloring problem

    // flag=0: non-specified; flag=1: upperNetwork; flag=2: lowerNetwork
    char* leftFlag = new char[size]();
    char* rightFlag = new char[size]();
    int rightPointer = size - 1;
    int leftPointer;
    while (rightFlag[rightPointer] == 0)
    {
        rightFlag[rightPointer] = 2;
        leftPointer = permuIndices[rightPointer];
        leftFlag[leftPointer] = 2;
        if (odd && leftPointer == size - 1)
            break;
        leftPointer = leftPointer & 1 ? leftPointer - 1 : leftPointer + 1;
        leftFlag[leftPointer] = 1;
        rightPointer = invPermuIndices[leftPointer];
        rightFlag[rightPointer] = 1;
        rightPointer = rightPointer & 1 ? rightPointer - 1 : rightPointer + 1;
    }
    for (int i = 0; i < size - 1; i++)
    {
        rightPointer = i;
        while (rightFlag[rightPointer] == 0)
        {
            rightFlag[rightPointer] = 2;
            leftPointer = permuIndices[rightPointer];
            leftFlag[leftPointer] = 2;
            leftPointer = leftPointer & 1 ? leftPointer - 1 : leftPointer + 1;
            leftFlag[leftPointer] = 1;
            rightPointer = invPermuIndices[leftPointer];

            rightFlag[rightPointer] = 1;
            rightPointer = rightPointer & 1 ? rightPointer - 1 : rightPointer + 1;
        }
    }
    delete[] invPermuIndices;

    // Determine bits on left gates
    int halfSize = size / 2;
    for (int i = 0; i < halfSize; i++)
        bits[i] = leftFlag[2 * i] == 2;

    int upperIndex = halfSize;
    int upperGateSize = getGateSize(halfSize);
    int lowerIndex = upperIndex + upperGateSize;
    int rightGateIndex = lowerIndex + (odd ? getGateSize(halfSize + 1) : upperGateSize);
    // Determine bits on right gates
    for (int i = 0; i < halfSize - 1; i++)
        bits[rightGateIndex + i] = rightFlag[2 * i] == 2;
    if (odd)
        bits[rightGateIndex + halfSize - 1] = rightFlag[size - 2] == 1;

    delete[] leftFlag;
    delete[] rightFlag;

    // Compute upper network
    int* upperIndices = new int[halfSize];
    for (int i = 0; i < halfSize - 1 + odd; i++)
        upperIndices[i] = permuIndices[2 * i + bits[rightGateIndex + i]] / 2;
    if (!odd)
        upperIndices[halfSize - 1] = permuIndices[size - 2] / 2;
    permutationToBits(upperIndices, halfSize, bits + upperIndex);
    delete[] upperIndices;

    // Compute lower network
    int lowerSize = halfSize + odd;
    int* lowerIndices = new int[lowerSize];
    for (int i = 0; i < halfSize - 1 + odd; i++)
        lowerIndices[i] = permuIndices[2 * i + 1 - bits[rightGateIndex + i]] / 2;
    if (odd)
        lowerIndices[halfSize] = permuIndices[size - 1] / 2;
    else
        lowerIndices[halfSize - 1] = permuIndices[2 * halfSize - 1] / 2;
    permutationToBits(lowerIndices, lowerSize, bits + lowerIndex);
    delete[] lowerIndices;

}

struct GateBlinder
{
    uint16_t upper;
    uint16_t lower;
};

// Inputs of the gate: v0=x0-r1, v1=x1-r2
// Outputs of the gate: if bit==1 then v0=x1-r3, v1=x0-r4; otherwise v0=x0-r3, v1=x1-r4
// m0=r1-r3, m1=r2-r4
void evaluateGate(uint16_t & v0, uint16_t & v1, GateBlinder blinder, bool bit)
{
    if (bit)
    {
        uint16_t temp = v1 + blinder.upper;
        v1 = v0 + blinder.lower;
        v0 = temp;
    }
    else
    {
        v0 += blinder.upper;
        v1 += blinder.lower;
    }
}


// If you want to apply the original exchange operation, set blinders to be 0;
void evaluateNetwork(uint16_t * values, int size, bool* bits, GateBlinder* blinders)
{
    if (size == 2)
        evaluateGate(values[0], values[1], blinders[0], bits[0]);
    if (size <= 2)
        return;

    int odd = size & 1;
    int halfSize = size / 2;

    // Compute left gates
    for (int i = 0; i < halfSize; i++)
        evaluateGate(values[2 * i], values[2 * i + 1], blinders[i], bits[i]);
    bits += halfSize;
    blinders += halfSize;

    // Compute upper subnetwork
    uint16_t* upperValues = new uint16_t[halfSize];
    for (int i = 0; i < halfSize; i++)
        upperValues[i] = values[i * 2];
    evaluateNetwork(upperValues, halfSize, bits, blinders);
    int upperGateSize = getGateSize(halfSize);
    bits += upperGateSize;
    blinders += upperGateSize;

    // Compute lower subnetwork
    int lowerSize = halfSize + odd;
    uint16_t* lowerValues = new uint16_t[lowerSize];
    for (int i = 0; i < halfSize; i++)
        lowerValues[i] = values[i * 2 + 1];
    if (odd) // the last element
        lowerValues[lowerSize - 1] = values[size - 1];
    evaluateNetwork(lowerValues, lowerSize, bits, blinders);
    int lowerGateSize = odd ? getGateSize(lowerSize) : upperGateSize;
    bits += lowerGateSize;
    blinders += lowerGateSize;

    // Deal with outputs of subnetworks
    for (int i = 0; i < halfSize; i++)
    {
        values[2 * i] = upperValues[i];
        values[2 * i + 1] = lowerValues[i];
    }
    if (odd) // the last element
        values[size - 1] = lowerValues[lowerSize - 1];

    // Compute right gates
    for (int i = 0; i < halfSize - 1 + odd; i++)
        evaluateGate(values[2 * i], values[2 * i + 1], blinders[i], bits[i]);

    delete[] upperValues;
    delete[] lowerValues;
}

struct Label
{
    uint16_t input1;
    uint16_t input2;
    uint16_t output1;
    uint16_t output2;
};

void writeGateLabels(uint16_t* inputLabel, int size, Label* gateLabels)
{
    if (size == 2)
    {
        gateLabels[0].input1 = inputLabel[0];
        gateLabels[0].input2 = inputLabel[1];
        gateLabels[0].output1 = gRNG.nextUInt16();
        gateLabels[0].output2 = gRNG.nextUInt16();
        inputLabel[0] = gateLabels[0].output1;
        inputLabel[1] = gateLabels[0].output2;
    }

    if (size <= 2)
        return;

    int odd = size & 1;
    int halfSize = size / 2;

    // Compute left gates
    for (int i = 0; i < halfSize; i++)
    {
        gateLabels[i].input1 = inputLabel[2 * i];
        gateLabels[i].input2 = inputLabel[2 * i + 1];
        gateLabels[i].output1 = gRNG.nextUInt16();
        gateLabels[i].output2 = gRNG.nextUInt16();
        inputLabel[2 * i] = gateLabels[i].output1;
        inputLabel[2 * i + 1] = gateLabels[i].output2;
    }
    gateLabels += halfSize;

    // Compute upper subnetwork
    uint16_t* upperInputs = new uint16_t[halfSize];
    for (int i = 0; i < halfSize; i++)
        upperInputs[i] = inputLabel[2 * i];
    writeGateLabels(upperInputs, halfSize, gateLabels);
    int upperGateSize = getGateSize(halfSize);
    gateLabels += upperGateSize;

    // Compute lower subnetwork
    int lowerSize = halfSize + odd;
    uint16_t* lowerInputs = new uint16_t[lowerSize];
    for (int i = 0; i < halfSize; i++)
        lowerInputs[i] = inputLabel[2 * i + 1];
    if (odd) // the last element
        lowerInputs[lowerSize - 1] = inputLabel[size - 1];
    writeGateLabels(lowerInputs, lowerSize, gateLabels);
    int lowerGateSize = odd ? getGateSize(lowerSize) : upperGateSize;
    gateLabels += lowerGateSize;

    // Deal with outputs of subnetworks
    for (int i = 0; i < halfSize; i++)
    {
        inputLabel[2 * i] = upperInputs[i];
        inputLabel[2 * i + 1] = lowerInputs[i];
    }
    if (odd) // the last element
        inputLabel[size - 1] = lowerInputs[lowerSize - 1];

    // Compute right gates
    for (int i = 0; i < halfSize - 1 + odd; i++)
    {
        gateLabels[i].input1 = inputLabel[2 * i];
        gateLabels[i].input2 = inputLabel[2 * i + 1];
        gateLabels[i].output1 = gRNG.nextUInt16();
        gateLabels[i].output2 = gRNG.nextUInt16();
        inputLabel[2 * i] = gateLabels[i].output1;
        inputLabel[2 * i + 1] = gateLabels[i].output2;
    }
    delete[] upperInputs;
    delete[] lowerInputs;
}


// oblivious permutation
// Alice receives Z1, while Bob computes Z2
void oblivPermutation(int* permutedIndices, uint16_t* values, int N, uint16_t* Z1, uint16_t* Z2)
{
    int gateSize = getGateSize(N);
    bool* bits = new bool[gateSize];

    // Alice generates selection bits
    permutationToBits(permutedIndices, N, bits);

    // Bob generates blinded inputs
    Label* gateLabels = new Label[gateSize];
    for (int i = 0; i < N; i++) 
        shareAnnot(values[i], Z1[i], Z2[i]); // Sends blinded inputs to Alice
    // Bob locally randomly writes labels for each gate
    writeGateLabels(Z2, N, gateLabels);

    // OT; Alice sets gateBlinders
    GateBlinder* gateBlinders = new GateBlinder[gateSize];
    uint16_t msg0[2], msg1[2], msgrecv[2];
    for (int i = 0; i < gateSize; i++)
    {
        Label label = gateLabels[i];
        msg0[0] = label.input1 - label.output1;
        msg0[1] = label.input2 - label.output2;
        msg1[0] = label.input2 - label.output1;
        msg1[1] = label.input1 - label.output2;
        OT::OT((uint8_t*)msg0, (uint8_t*)msg1, 16/8*2, bits[i], (uint8_t*)msgrecv);
        gateBlinders[i].upper = msgrecv[0];
        gateBlinders[i].lower = msgrecv[1];
    }

    // Alice evalutes the network
    evaluateNetwork(Z1, N, bits, gateBlinders);

    delete[] bits;
    delete[] gateLabels;
    delete[] gateBlinders;

}

void OEP(int* indices, int M, int N, uint16_t* values, uint16_t* Z1, uint16_t* Z2)
{
    int origM = M;
    if (N > M)
        M = N;
    uint16_t* extendedValues = new uint16_t[M];
    for (int i = 0; i < origM; i++)
        extendedValues[i] = values[i];
    // remained values are not assgined, which are not important
    int* indicesCount = new int[M]();
    for (int i = 0; i < N; i++)
        indicesCount[indices[i]]++;
    int* firstPermu = new int[M];
    int dummyIndex = 0, fPIndex = 0;
    // We call those index with indicesCount[index]==0 as dummy index
    for(int i=0;i<M;i++)
    {
        if (indicesCount[i] > 0)
        {
            firstPermu[fPIndex++] = i;
            for (int j = 0; j < indicesCount[i] - 1; j++)
            {
                while (indicesCount[dummyIndex] > 0)
                    dummyIndex++;
                firstPermu[fPIndex++] = dummyIndex++;
            }
        }
    }
    while (fPIndex < M)
    {
        while (indicesCount[dummyIndex] > 0)
            dummyIndex++;
        firstPermu[fPIndex++] = dummyIndex++;
    }

    uint16_t* Z11 = new uint16_t[M];
    uint16_t* Z12 = new uint16_t[M];


    oblivPermutation(firstPermu, extendedValues, M, Z11, Z12);

    // Replication network
    for (int i = 1; i < N; i++)
    {
        // Alice determines the bit
        bool bit = indicesCount[firstPermu[i]] == 0; // dummy index

        // Bob computes the messages 
        uint16_t r = gRNG.nextUInt16();
        uint16_t msg0 = Z12[i] - r;
        uint16_t msg1 = Z12[i - 1] - r;
        uint16_t msgb = 0;

        // OT
        OT::OT((uint8_t*)&msg0, (uint8_t*)&msg1, 16/8, bit, (uint8_t*)&msgb);

        // Update the shares
        Z11[i] = Z11[i - bit] + msgb;
        Z12[i] = r;
    }

    int* pointers = new int[M];
    int sum = 0;
    for (int i = 0; i < M; i++)
    {
        pointers[i] = sum;
        sum += indicesCount[i];
    }
    int* totalMap = new int[N];
    for (int i = 0; i < N; i++)
        totalMap[i] = firstPermu[pointers[indices[i]]++];
    int* invFirstPermu = new int[M];
    for (int i = 0; i < M; i++)
        invFirstPermu[firstPermu[i]] = i;
    int* secondPermu = new int[N];
    for (int i = 0; i < N; i++)
        secondPermu[i] = invFirstPermu[totalMap[i]];

    uint16_t* Z21 = new uint16_t[N];
    uint16_t* Z22 = new uint16_t[N];
    oblivPermutation(secondPermu, Z12, N, Z21, Z22);
    for (int i = 0; i < N; i++)
    {
        Z1[i] = Z11[secondPermu[i]] + Z21[i];
        Z2[i] = Z22[i];
    }

    delete[] indicesCount;
    delete[] firstPermu;
    delete[] invFirstPermu;
    delete[] extendedValues;
    delete[] secondPermu;
    delete[] totalMap;
    delete[] pointers;
    delete[] Z11;
    delete[] Z12;
    delete[] Z21;
    delete[] Z22;

}