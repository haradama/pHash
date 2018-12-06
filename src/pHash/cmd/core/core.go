package core

import (
	"runtime"
	"strings"

	"github.com/OneOfOne/xxhash"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/seq"
)

type Plasmids struct {
	Plasmid map[string][]uint64
}

func revComp(seq string) string {

	ambiguousDnaComplement := strings.NewReplacer(
		"A", "T",
		"C", "G",
		"G", "C",
		"T", "A",
		"M", "K",
		"R", "Y",
		"Y", "R",
		"K", "M",
		"V", "B",
		"H", "D",
		"D", "H",
		"B", "V")

	runes := []rune(seq)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}

	seqRC := ambiguousDnaComplement.Replace(string(runes))

	return seqRC
}

func MinHash(s *seq.Sequence, k int, NGoRoutines int, sketchSize int) *[]uint64 {
	read := (*s).Slice()
	kmerNum := read.Len() - (k - 1)
	kmerMap := make(map[string]struct{})
	for i := 0; i < kmerNum; i++ {
		seq := read.Slice(i, i+k).(alphabet.Letters).String()
		if strings.Index(seq, "N") == -1 {
			reverseComplement := revComp(seq)

			h1 := xxhash.ChecksumString32S(seq, 0)
			h2 := xxhash.ChecksumString32S(reverseComplement, 0)

			var canonicalKmer string
			if h1 > h2 {
				canonicalKmer = seq
			} else {
				canonicalKmer = reverseComplement
			}
			kmerMap[canonicalKmer] = struct{}{}
		}
	}

	kmerList := [][]byte{}
	for key := range kmerMap {
		kmerList = append(kmerList, []byte(key))
	}

	c := make(chan int, NGoRoutines)
	runtime.GOMAXPROCS(NGoRoutines)
	chunk := sketchSize / NGoRoutines

	minHashValues := make([]uint64, sketchSize)

	for i := 0; i < NGoRoutines; i++ {
		go func(start int) {
			end := start + chunk

			if end > sketchSize {
				end = sketchSize
			}

			for j := uint64(start); j < uint64(end); j++ {
				minValue := uint64((1 << 64) - 1)
				for _, kmer := range kmerList {
					xxhv := xxhash.Checksum64S(kmer, j)
					if minValue > xxhv {
						minValue = xxhv
					}
				}
				minHashValues[j] = minValue
			}
			c <- 1
		}(i * chunk)
	}

	for i := 0; i < NGoRoutines; i++ {
		<-c
	}

	return &minHashValues
}
