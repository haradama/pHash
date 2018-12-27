package cmd

import (
	"bufio"
	"log"
	"os"
	"strings"
	"sync"

	"github.com/OneOfOne/xxhash"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/seq"
	"github.com/ugorji/go/codec"
	pb "gopkg.in/cheggaaa/pb.v1"
)

func makeBar(inFile *string) *pb.ProgressBar {
	f, err := os.Open(*inFile)
	if err != nil {
		os.Exit(1)
	}
	header := byte('>')
	var counter int

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		text := scanner.Text()
		if len(text) > 0 && text[0] == header {
			counter++
		}
	}

	bar := pb.New(counter)
	bar.ShowPercent = false
	bar.ShowSpeed = false
	bar.ShowTimeLeft = false
	bar.ShowBar = false
	bar.ShowCounters = false

	return bar
}

func messagePackEncoding(plasmids *Plasmids) []byte {
	buf := make([]byte, 0, 64)
	err := codec.NewEncoderBytes(&buf, &mh).Encode(*plasmids)
	if err != nil {
		log.Printf("error encoding %v to MessagePack: %v", *plasmids, err)
	}
	return buf
}

func messagePackDecoding(buf *[]byte) Plasmids {
	var plasmids Plasmids
	err := codec.NewDecoderBytes(*buf, &mh).Decode(&plasmids)
	if err != nil {
		log.Printf("error decoding %v to MessagePack: %v", buf, err)
	}
	return plasmids
}

func getKmerList(s *seq.Sequence, k *int, ambiguousDnaComplement *strings.Replacer) *[][]byte {
	mutex := &sync.Mutex{}
	wg := sync.WaitGroup{}

	read := (*s).Slice()
	kmerNum := read.Len() - (*k - 1)
	kmerMap := make(map[string]struct{})

	for i := 0; i < kmerNum; i++ {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			seq := read.Slice(i, i+*k).(alphabet.Letters).String()
			if strings.Index(seq, "N") == -1 {
				reverseComplement := ambiguousDnaComplement.Replace(rev(&seq))

				h1 := xxhash.ChecksumString32S(seq, 0)
				h2 := xxhash.ChecksumString32S(reverseComplement, 0)

				var canonicalKmer string
				if h1 > h2 {
					canonicalKmer = seq
				} else {
					canonicalKmer = reverseComplement
				}
				mutex.Lock()
				kmerMap[canonicalKmer] = struct{}{}
				mutex.Unlock()
			}
		}(i)
	}
	wg.Wait()

	kmerList := [][]byte{}
	for key := range kmerMap {
		kmerList = append(kmerList, []byte(key))
	}

	return &kmerList
}
