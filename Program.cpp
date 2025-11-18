#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

using namespace std;

// Parameters
const int L = 12;         // Segment length
const int K = 3;          // k-mer size
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
vector<Segment> SegmentSequence(const string &dna);
string ComputeRepresentativeWord(const string &segment);
vector<double> ComputePositionPValues(const string &segment);
double BinomialPValue(int n, int k, double p);
double BinomialCoefficient(int n, int k);
double CombinePValuesFisher(const vector<double> &pValues);
vector<Segment> MergeSameWordSegments(const vector<Segment> &segments);
vector<Segment> MergeNoiseSegments(const vector<Segment> &segments);

int main()
{
    // FIX: Removed '+' signs to allow implicit string literal concatenation
    string dna =
        "GTGACGGTGTAG"  // strong repeat GTG
        "ACGTTAGGACTA"  // weak noise
        "GTGACGGTGTAG"; // strong repeat GTG again

    cout << "Loaded sequence length: " << dna.length() << endl;
    cout << "DNA Sequence Preview: " << dna << endl;

    // 1. Segment sequence
    vector<Segment> segments = SegmentSequence(dna);

    // 2. Merge same-word segments
    segments = MergeSameWordSegments(segments);

    // 3. Merge weak/noisy segments
    segments = MergeNoiseSegments(segments);

    // Output results
    cout << "\nDetected Hidden Repeat Segments:\n";
    for (auto &s : segments)
    {
        cout << "Start: " << s.startIndex
             << ", Len: " << s.length
             << ", Word: " << s.representativeWord
             << ", Score(P): " << s.combinedPValue
             << endl;
    }

    return 0;
}

// Segment sequence into chunks and compute words + p-values
vector<Segment> SegmentSequence(const string &dna)
{
    vector<Segment> segments;

    for (size_t i = 0; i < dna.size(); i += L)
    {
        string subseq = dna.substr(i, min((size_t)L, dna.size() - i));
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
string ComputeRepresentativeWord(const string &segment)
{
    string word = "";
    int numKmers = segment.size() / K;

    for (int pos = 0; pos < K; ++pos)
    {
        map<char, int> counts{{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};

        for (int i = 0; i < numKmers; ++i)
        {
            int idx = i * K + pos;
            if (idx >= segment.size())
                continue;
            counts[segment[idx]]++;
        }
        // Find max
        char maxNuc = 'A';
        int maxCount = 0;
        for (auto &kv : counts)
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
vector<double> ComputePositionPValues(const string &segment)
{
    vector<double> pValues(K);
    int numKmers = segment.size() / K;

    for (int pos = 0; pos < K; ++pos)
    {
        map<char, int> counts{{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};

        for (int i = 0; i < numKmers; ++i)
        {
            int idx = i * K + pos;
            if (idx >= segment.size())
                continue;
            counts[segment[idx]]++;
        }
        int maxCount = 0;
        for (auto &kv : counts)
            maxCount = max(maxCount, kv.second);
        pValues[pos] = BinomialPValue(numKmers, maxCount, 0.25);
    }
    return pValues;
}

// Binomial P(X>=k)
double BinomialPValue(int n, int k, double p)
{
    double sum = 0.0;
    for (int i = k; i <= n; ++i)
        sum += BinomialCoefficient(n, i) * pow(p, i) * pow(1.0 - p, n - i);
    return sum;
}

// Binomial coefficient nCk
double BinomialCoefficient(int n, int k)
{
    if (k < 0 || k > n)
        return 0;
    if (k == 0 || k == n)
        return 1;
    if (k > n / 2)
        k = n - k;

    double res = 1.0;
    for (int i = 1; i <= k; ++i)
        res = res * (n - i + 1) / i;
    return res;
}

// Combine p-values using Fisher's method
double CombinePValuesFisher(const vector<double> &pValues)
{
    double X = 0.0;
    for (double p : pValues)
    {
        // Avoid log(0) if p is extremely small
        if (p > 1e-300)
            X += -2.0 * log(p);
        else
            X += -2.0 * log(1e-300);
    }
    // Note: This simplified return acts as a score, not a direct Chi-sq p-value
    return exp(-0.5 * X);
}

// Merge consecutive segments with same word
vector<Segment> MergeSameWordSegments(const vector<Segment> &segments)
{
    vector<Segment> merged;
    if (segments.empty())
        return merged;

    Segment current = segments[0];

    for (size_t i = 1; i < segments.size(); ++i)
    {
        if (segments[i].representativeWord == current.representativeWord)
        {
            current.sequence += segments[i].sequence;
            current.length += segments[i].length;
            // Recompute stats for the merged segment
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
vector<Segment> MergeNoiseSegments(const vector<Segment> &segments)
{
    vector<Segment> result;
    if (segments.empty())
        return result;

    // We need a copy of segments because we might modify flow,
    // but here we are just iterating.
    // Since we are "skipping" indices, a while loop is appropriate.

    result.push_back(segments[0]); // Start with the first

    size_t i = 1;
    while (i < segments.size())
    {
        // We look at the *last pushed* segment (left) and the current one (middle)
        // We need to check if there is a 'right' segment to bridge to.
        if (i < segments.size() - 1)
        {
            Segment &left = result.back();
            const Segment &middle = segments[i];
            const Segment &right = segments[i + 1];

            bool leftStrong = left.combinedPValue < tau2;
            bool rightStrong = right.combinedPValue < tau2;
            bool middleWeak = middle.combinedPValue > tau1; // High p-value = weak signal

            // Check if middle is noise between two strong signals
            // Note: Original logic checked words weren't equal to middle,
            // implying middle is an interruption.
            if (middle.representativeWord != left.representativeWord &&
                middle.representativeWord != right.representativeWord &&
                leftStrong && rightStrong && middleWeak)
            {
                // Merge Middle and Right into Left
                left.sequence += middle.sequence + right.sequence;
                left.length = left.sequence.size();
                left.positionPValues = ComputePositionPValues(left.sequence);
                left.combinedPValue = CombinePValuesFisher(left.positionPValues);

                i += 2; // Skip middle and right
                continue;
            }
        }
        // Otherwise just add the current segment
        result.push_back(segments[i]);
        i++;
    }
    return result;
}