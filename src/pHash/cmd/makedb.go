package cmd

import (
	"fmt"
	"io"
	"os"
	"strings"
	"sync"

	"github.com/OneOfOne/xxhash"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/spf13/cobra"
)

func init() {
	RootCmd.AddCommand(makedbCmd)
	makedbCmd.Flags().StringVarP(&o.optIn, "in", "i", "default", "Input FASTA file")
	makedbCmd.Flags().StringVarP(&o.optBuildOut, "out", "o", "reference.phash", "Database")
	makedbCmd.Flags().IntVarP(&o.optKmer, "kmer", "k", 16, "Length of k-mer")
	makedbCmd.Flags().IntVarP(&o.optSketch, "sketch", "s", 1024, "Sketch size")
}

var makedbCmd = &cobra.Command{
	Use:   "makedb",
	Short: "Builder of plasmid database",
	Long:  "Builder of plasmid database using MinHash",
	Run: func(cmd *cobra.Command, args []string) {

		if o.optIn == "" || o.optBuildOut == "" {
			fmt.Println("--in and --out are required")
			cmd.Help()
			os.Exit(0)
		}

		inFile := o.optIn
		outFile := o.optBuildOut
		k := o.optKmer
		sketchSize := o.optSketch

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

		var wg sync.WaitGroup
		mutex := new(sync.Mutex)

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

		plasmidsMap := map[string][]uint64{}

		for i := 0; i < 1000; i++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				for {
					mutex.Lock()
					s, err := in.Read()
					mutex.Unlock()

					if err != nil {
						if err != io.EOF {
							fmt.Println(err)
							os.Exit(1)
						}
						break
					}

					read := s.Slice()
					kmerNum := read.Len() - (k - 1)
					kmerMap := make(map[string]struct{})

					for i := 0; i < kmerNum; i++ {
						seq := read.Slice(i, i+k).(alphabet.Letters).String()
						if strings.Index(seq, "N") == -1 {
							reverseComplement := ambiguousDnaComplement.Replace(rev(&seq))
							h1 := getCanonicalKmer(seq)
							h2 := getCanonicalKmer(reverseComplement)

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

					minHashValues := make([]uint64, sketchSize)

					for i := uint64(0); i < uint64(sketchSize); i++ {
						minValue := uint64((1 << 64) - 1)
						for _, kmer := range kmerList {
							xxhv := xxhash.Checksum64S(kmer, i)
							if minValue > xxhv {
								minValue = xxhv
							}
						}
						minHashValues[i] = minValue
					}

					mutex.Lock()
					plasmidsMap[s.Name()] = minHashValues
					mutex.Unlock()
				}
			}()
		}
		wg.Wait()

		plasmids := Plasmids{
			SketchSize: sketchSize,
			Kmer:       k,
			Plasmid:    plasmidsMap,
		}

		file, err := os.Create(outFile)
		if err != nil {
			if err != io.EOF {
				fmt.Println(err)
				os.Exit(1)
			}
		}
		defer file.Close()

		file.Write(messagePackEncoding(&plasmids))
	},
}
