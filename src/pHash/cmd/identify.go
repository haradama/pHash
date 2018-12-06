package cmd

import (
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"math"
	"os"

	"./core"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/spf13/cobra"
	"github.com/ugorji/go/codec"
)

func init() {
	RootCmd.AddCommand(identifyCmd)
	identifyCmd.Flags().StringVarP(&o.optIn, "in", "i", "", "Input FASTA file")
	identifyCmd.Flags().StringVarP(&o.optIdentifyOut, "out", "o", "plasmid.fna", "Output FASTA file")
	identifyCmd.Flags().StringVarP(&o.optDB, "db", "d", "plasmidDB.phash", "Database")
	identifyCmd.Flags().IntVarP(&o.optThreshold, "threshold", "t", 10, "Threshold of probability")
	identifyCmd.Flags().IntVarP(&o.optParalell, "paralell", "p", 4, "Number of parallel processing")
	identifyCmd.Flags().IntVarP(&o.optKmer, "kmer", "k", 17, "Length of k-mer")
	identifyCmd.Flags().IntVarP(&o.optSketch, "sketch", "s", 1024, "Sketch size")
}

var identifyCmd = &cobra.Command{
	Use:   "identify",
	Short: "Identifier of plasmid",
	Long:  "Identifier of plasmid using database",
	Run: func(cmd *cobra.Command, args []string) {
		inFile := o.optIn
		db := o.optDB
		outFile := o.optIdentifyOut
		paralell := o.optParalell
		k := o.optKmer
		sketchSize := o.optSketch
		threshold := float64(o.optThreshold) * 0.01

		var in *fasta.Reader
		if inFile == "" {
			in = fasta.NewReader(os.Stdin, linear.NewSeq("", nil, alphabet.DNA))
		} else if f, err := os.Open(inFile); err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v.", err)
			os.Exit(1)
		} else {
			in = fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNA))
			defer f.Close()
		}

		binary, err := ioutil.ReadFile(db)
		if err != nil {
			fmt.Println("File reading error", err)
			return
		}

		plasmidsMap := messagePackDecoding(binary).Plasmid

		file, _ := os.Create(outFile)
		w := fasta.NewWriter(file, 60)
		plasmidSeq := linear.NewSeq("", nil, alphabet.DNA)

		for {
			s, err := in.Read()
			if err != nil {
				if err != io.EOF {
					fmt.Println(err)
					os.Exit(1)
				}
				break
			}
			minHashValue := core.MinHash(&s, k, paralell, sketchSize)

			var (
				bestHitKey   string
				bestHitValue float64 = 0
			)
			similarityMap := make(map[string]float64)

			for refKey, refValue := range plasmidsMap {
				similarity := calcSimilarity(*minHashValue, refValue, 1024)
				similarityMap[refKey] = similarity

				if similarity > bestHitValue {
					bestHitKey = refKey
					bestHitValue = similarity
				}
			}
			if bestHitValue >= threshold {
				plasmidSeq.ID = s.Name()
				plasmidSeq.Desc = fmt.Sprintf("Similar to %s (%f)", bestHitKey, bestHitValue)
				_s := s.Slice()
				plasmidSeq.Seq = []alphabet.Letter(_s.(alphabet.Letters).String())
				fmt.Println(s.Name(), bestHitKey, bestHitValue)
				_, err := w.Write(plasmidSeq)
				if err != nil {
					log.Printf("Failed to write: %s", err)
				}
			}
		}
	},
}

func messagePackDecoding(buf []byte) Plasmids {
	var plasmids Plasmids
	err := codec.NewDecoderBytes(buf, &mh).Decode(&plasmids)
	if err != nil {
		log.Printf("error decoding %v to MessagePack: %v", buf, err)
	}
	return plasmids
}

func calcSimilarity(array1 []uint64, array2 []uint64, size float64) float64 {
	var match float64
	for i := 0; i < int(size); i++ {
		if array1[i] == array2[i] {
			match++
		}
	}
	var bitNum float64 = 64
	similarity := (match / size) - math.Pow(2, -bitNum)/(1-math.Pow(2, -bitNum))

	return similarity
}
