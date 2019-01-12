package cmd

import (
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"runtime"
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
	makedbCmd.Flags().StringVarP(&o.optMetadata, "meta", "m", "", "Input FASTA file")
	makedbCmd.Flags().StringVarP(&o.optBuildOut, "out", "o", "reference.phash", "Database")
	makedbCmd.Flags().IntVarP(&o.optKmer, "kmer", "k", 16, "Length of k-mer")
	makedbCmd.Flags().IntVarP(&o.optSketch, "sketch", "s", 512, "Sketch size")
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
		sketchSize := uint64(o.optSketch)
		metadata := o.optMetadata

		cpus := runtime.NumCPU()
		runtime.GOMAXPROCS(cpus)

		phylumMap := map[string]string{}
		if len(metadata) > 0 {
			csvfile, err := os.Open(metadata)
			if err != nil {
				if err != io.EOF {
					fmt.Println(err)
					os.Exit(1)
				}
				return
			}
			defer csvfile.Close()

			csvreader := csv.NewReader(csvfile)

			if metadata != "" {
				for {
					record, _ := csvreader.Read()
					if err != nil {
						if err != io.EOF {
							fmt.Println(err)
							os.Exit(1)
						}
					}
					if len(record) == 0 {
						break
					}
					phylumMap[record[0]] = record[1]
				}
			}
		}

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

		plasmidsRecords := []PlasmidRecord{}

		for i := 0; i < 1024; i++ {
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
					kmerMap := make(map[string]struct{}, kmerNum)

					for i := 0; i < kmerNum; i++ {
						seq := read.Slice(i, i+k).(alphabet.Letters).String()
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
							kmerMap[canonicalKmer] = struct{}{}
						}
					}

					kmerList := make([][]byte, len(kmerMap))
					for key := range kmerMap {
						kmerList = append(kmerList, []byte(key))
					}

					minHashValues := make([]uint64, sketchSize)

					for i := uint64(0); i < sketchSize; i++ {
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
					phylum := "---"
					if value, ok := phylumMap[s.Name()]; ok {
						phylum = value
					}
					plasmidsRecords = append(plasmidsRecords, PlasmidRecord{AccID: s.Name(), Phylum: phylum, PlasmidMinHashValue: minHashValues})
					mutex.Unlock()
				}
			}()
		}
		wg.Wait()

		plasmids := Plasmids{
			SketchSize: sketchSize,
			Kmer:       k,
			Plasmid:    plasmidsRecords,
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
