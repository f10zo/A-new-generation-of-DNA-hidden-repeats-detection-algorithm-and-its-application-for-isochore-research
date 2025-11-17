#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

using namespace std;

// Parameters
const int L = 12;        // Segment length
const int K = 3;         // k-mer size
const double tau1 = 0.1;  // High p-value threshold
const double tau2 = 0.01; // Low p-value threshold
const double alpha = 0.05;

struct Segment
{
    string sequence;
    string representativeWord;
    vector<double> positionPValues;
    double combinedPValue;
    int startIndex;
    int length;
};

// Function declarations
vector<Segment> SegmentSequence(const string& dna);
string ComputeRepresentativeWord(const string& segment);
vector<double> ComputePositionPValues(const string& segment);
double BinomialPValue(int n, int k, double p);
double BinomialCoefficient(int n, int k);
double CombinePValuesFisher(const vector<double>& pValues);
vector<Segment> MergeSameWordSegments(const vector<Segment>& segments);
vector<Segment> MergeNoiseSegments(const vector<Segment>& segments);

int main()
{
    string dna =
        "GTGACGGTGTAG" +   // strong repeat GTG
        "ACGTTAGGACTA" +   // weak noise
        "GTGACGGTGTAG";    // strong repeat GTG again → should merge

    cout << "Loaded sequence length: " << dna.length() << endl;
    cout << "DNA Sequence Preview: " << dna << endl;

    // Segment sequence
    vector<Segment> segments = SegmentSequence(dna);

    // Merge same-word segments
    segments = MergeSameWordSegments(segments);

    // Merge weak/noisy segments
    segments = MergeNoiseSegments(segments);

    // Output results
    cout << "\nDetected Hidden Repeat Segments:\n";
    for (auto & s : segments)
    {
        cout << "Start:" << s.startIndex
             << ", Len:" << s.length
             << ", Word:" << s.representativeWord
             << ", P:" << s.combinedPValue
             << endl;
    }

    return 0;
}

// Segment sequence into chunks and compute words + p-values
vector<Segment> SegmentSequence(const string& dna)
{
    vector<Segment> segments;

    for (size_t i = 0; i < dna.size(); i += L)
    {
        string subseq = dna.substr(i, min(L, (int)(dna.size() - i)));
        Segment seg;
        seg.sequence = subseq;
        seg.representativeWord = ComputeRepresentativeWord(subseq);
        seg.positionPValues = ComputePositionPValues(subseq);
        seg.startIndex = i;
        seg.length = subseq.size();
        seg.combinedPValue = CombinePValuesFisher(seg.positionPValues);
        segments.push_back(seg);
    }

    return segments;
}

// Compute representative k-mer word
string ComputeRepresentativeWord(const string& segment)
{
    string word = "";
    int numKmers = segment.size() / K;

    for (int pos = 0; pos < K; ++pos)
    {
        map<char, int> counts{ { 'A',0},{ 'C',0},{ 'G',0},{ 'T',0} }
        ;
        for (int i = 0; i < numKmers; ++i)
        {
            int idx = i * K + pos;
            if (idx >= segment.size()) continue;
            counts[segment[idx]]++;
        }
        // Find max
        char maxNuc = 'A';
        int maxCount = 0;
        for (auto & kv : counts)
        {
            if (kv.second > maxCount)
            {
                maxCount = kv.second;
                maxNuc = kv.first;
            }
        }
        word += maxNuc;
    }

    return word;
}

// Compute p-values for each position
vector<double> ComputePositionPValues(const string& segment)
{
    vector<double> pValues(K);
    int numKmers = segment.size() / K;

    for (int pos = 0; pos < K; ++pos)
    {
        map<char, int> counts{ { 'A',0},{ 'C',0},{ 'G',0},{ 'T',0} }
        ;
        for (int i = 0; i < numKmers; ++i)
        {
            int idx = i * K + pos;
            if (idx >= segment.size()) continue;
            counts[segment[idx]]++;
        }
        int maxCount = 0;
        for (auto & kv : counts) maxCount = max(maxCount, kv.second);
        pValues[pos] = BinomialPValue(numKmers, maxCount, 0.25);
    }
    return pValues;
}

// Binomial P(X>=k)
double BinomialPValue(int n, int k, double p)
{
    double sum = 0.0;
    for (int i = k; i <= n; ++i)
        sum += BinomialCoefficient(n, i) * pow(p, i) * pow(1 - p, n - i);
    return sum;
}

// Binomial coefficient nCk
double BinomialCoefficient(int n, int k)
{
    double res = 1.0;
    for (int i = 1; i <= k; ++i)
        res *= (double)(n - i + 1) / i;
    return res;
}

// Combine p-values using Fisher's method
double CombinePValuesFisher(const vector<double>& pValues)
{
    double X = 0.0;
    for (double p : pValues) if (p > 0) X += -2 * log(p);
    return exp(-0.5 * X); // simplified survival function
}

// Merge consecutive segments with same word
vector<Segment> MergeSameWordSegments(const vector<Segment>& segments)
{
    vector<Segment> merged;
    Segment current;

    for (size_t i = 0; i < segments.size(); ++i)
    {
        if (i == 0)
        {
            current = segments[i];
        }
        else if (segments[i].representativeWord == current.representativeWord)
        {
            current.sequence += segments[i].sequence;
            current.length += segments[i].length;
            current.positionPValues = ComputePositionPValues(current.sequence);
            current.combinedPValue = CombinePValuesFisher(current.positionPValues);
        }
        else
        {
            merged.push_back(current);
            current = segments[i];
        }
    }
    merged.push_back(current);
    return merged;
}

// Merge noisy segments
vector<Segment> MergeNoiseSegments(const vector<Segment>& segments)
{
    vector<Segment> result;
    size_t i = 0;
    while (i < segments.size())
    {
        if (i > 0 && i < segments.size() - 1)
        {
            Segment & left = result.back();
            const Segment &middle = segments[i];
            const Segment &right = segments[i + 1];

            bool leftStrong = left.combinedPValue < tau2;
            bool rightStrong = right.combinedPValue < tau2;
            bool middleWeak = middle.combinedPValue > tau1;

            if (middle.representativeWord != left.representativeWord &&
                middle.representativeWord != right.representativeWord &&
                leftStrong && rightStrong && middleWeak)
            {

                left.sequence += middle.sequence + right.sequence;
                left.length = left.sequence.size();
                left.positionPValues = ComputePositionPValues(left.sequence);
                left.combinedPValue = CombinePValuesFisher(left.positionPValues);

                i += 2; // skip middle & right
                continue;
            }
        }
        result.push_back(segments[i]);
        i++;
    }
    return result;
}