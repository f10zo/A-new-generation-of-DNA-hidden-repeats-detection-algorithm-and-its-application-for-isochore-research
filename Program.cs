using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

class Program
{
    // Parameters
    const int L = 12;        // Segment length
    const int K = 3;         // k-mer size
    const double tau1 = 0.1;  // High p-value threshold
    const double tau2 = 0.01; // Low p-value threshold
    const double alpha = 0.05;

    static void Main()
    {
        // 🔹 Short artificial DNA sequence with hidden repeated regions
        string dna =
            "GTGACGGTGTAG" +   // strong repeat GTG
            "ACGTTAGGACTA" +   // weak noise
            "GTGACGGTGTAG";    // strong repeat GTG again → should merge

        Console.WriteLine($"Loaded sequence length: {dna.Length}");
        Console.WriteLine("DNA Sequence Preview:");
        Console.WriteLine(dna);

        // Step 2–4: Segment and compute words + p-values
        var segments = SegmentSequence(dna, L, K);

        // Step 6: Merge identical words
        segments = MergeSameWordSegments(segments);

        // Step 7: Merge weak/noisy segments
        segments = MergeNoiseSegments(segments, tau1, tau2, alpha);

        // Output results
        Console.WriteLine("\nDetected Hidden Repeat Segments:");
        foreach (var s in segments)
        {
            Console.WriteLine(
                $"Start:{s.StartIndex}, Len:{s.Length}, Word:{s.RepresentativeWord}, P:{s.CombinedPValue:F5}");
        }
    }

    class Segment
    {
        public string Sequence;
        public string RepresentativeWord;
        public double[] PositionPValues;
        public double CombinedPValue;
        public int StartIndex;
        public int Length;
    }

    static List<Segment> SegmentSequence(string dna, int L, int K)
    {
        var segments = new List<Segment>();

        for (int i = 0; i < dna.Length; i += L)
        {
            string subseq = dna.Substring(i, Math.Min(L, dna.Length - i));
            var seg = new Segment
            {
                Sequence = subseq,
                RepresentativeWord = ComputeRepresentativeWord(subseq, K),
                PositionPValues = ComputePositionPValues(subseq, K),
                StartIndex = i,
                Length = subseq.Length
            };
            seg.CombinedPValue = CombinePValuesFisher(seg.PositionPValues);
            segments.Add(seg);
        }
        return segments;
    }

    static string ComputeRepresentativeWord(string segment, int k)
    {
        var sb = new StringBuilder();
        int numKmers = segment.Length / k;

        for (int pos = 0; pos < k; pos++)
        {
            var counts = new Dictionary<char, int> {
                { 'A', 0 }, { 'C', 0 }, { 'G', 0 }, { 'T', 0 }
            };

            for (int i = 0; i < numKmers; i++)
            {
                int idx = i * k + pos;
                if (idx >= segment.Length) continue;
                char c = segment[idx];
                if (counts.ContainsKey(c)) counts[c]++;
            }

            char maxNuc = counts.OrderByDescending(x => x.Value).First().Key;
            sb.Append(maxNuc);
        }
        return sb.ToString();
    }

    static double[] ComputePositionPValues(string segment, int k)
    {
        int numKmers = segment.Length / k;
        var pValues = new double[k];

        for (int pos = 0; pos < k; pos++)
        {
            var counts = new Dictionary<char, int> {
                { 'A', 0 }, { 'C', 0 }, { 'G', 0 }, { 'T', 0 }
            };

            for (int i = 0; i < numKmers; i++)
            {
                int idx = i * k + pos;
                if (idx >= segment.Length) continue;
                char c = segment[idx];
                if (counts.ContainsKey(c)) counts[c]++;
            }

            int maxCount = counts.Values.Max();
            pValues[pos] = BinomialPValue(numKmers, maxCount, 0.25);
        }
        return pValues;
    }

    static double BinomialPValue(int n, int k, double p)
    {
        double sum = 0.0;
        for (int i = k; i <= n; i++)
        {
            sum += BinomialCoefficient(n, i) * Math.Pow(p, i) * Math.Pow(1 - p, n - i);
        }
        return sum;
    }

    static double BinomialCoefficient(int n, int k)
    {
        double res = 1;
        for (int i = 1; i <= k; i++)
            res *= (double)(n - i + 1) / i;
        return res;
    }

    static double CombinePValuesFisher(double[] pValues)
    {
        double X = -2 * pValues.Where(p => p > 0).Sum(p => Math.Log(p));
        return Math.Exp(-0.5 * X);
    }

    static List<Segment> MergeSameWordSegments(List<Segment> segments)
    {
        var merged = new List<Segment>();
        Segment current = null;

        foreach (var seg in segments)
        {
            if (current == null)
            {
                current = seg;
            }
            else if (current.RepresentativeWord == seg.RepresentativeWord)
            {
                current.Sequence += seg.Sequence;
                current.Length += seg.Length;
                current.PositionPValues = ComputePositionPValues(current.Sequence, K);
                current.CombinedPValue = CombinePValuesFisher(current.PositionPValues);
            }
            else
            {
                merged.Add(current);
                current = seg;
            }
        }
        if (current != null) merged.Add(current);
        return merged;
    }

    static List<Segment> MergeNoiseSegments(List<Segment> segments, double tau1, double tau2, double alpha)
    {
        var result = new List<Segment>();
        int i = 0;

        while (i < segments.Count)
        {
            if (i > 0 && i < segments.Count - 1)
            {
                var left = result.Last();
                var middle = segments[i];
                var right = segments[i + 1];

                bool leftStrong = left.CombinedPValue < tau2;
                bool rightStrong = right.CombinedPValue < tau2;
                bool middleWeak = middle.CombinedPValue > tau1;

                if (middle.RepresentativeWord != left.RepresentativeWord &&
                    middle.RepresentativeWord != right.RepresentativeWord &&
                    leftStrong && rightStrong && middleWeak)
                {
                    left.Sequence += middle.Sequence + right.Sequence;
                    left.Length = left.Sequence.Length;
                    left.PositionPValues = ComputePositionPValues(left.Sequence, K);
                    left.CombinedPValue = CombinePValuesFisher(left.PositionPValues);

                    i += 2; // Skip middle + right
                    continue;
                }
            }

            result.Add(segments[i]);
            i++;
        }

        return result;
    }
}
