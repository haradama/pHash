package cmd

import (
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"math"
	"os"
	"runtime"
	"strings"
	"sync"

	"github.com/OneOfOne/xxhash"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/spf13/cobra"
	pb "gopkg.in/cheggaaa/pb.v1"
)

func init() {
	RootCmd.AddCommand(identifyCmd)
	identifyCmd.Flags().StringVarP(&o.optIn, "in", "i", "", "Input FASTA file")
	identifyCmd.Flags().StringVarP(&o.optIdentifyOut, "out", "o", "", "Output FASTA file")
	identifyCmd.Flags().StringVarP(&o.optDB, "db", "d", "", "Database")
	identifyCmd.Flags().IntVarP(&o.optThreshold, "threshold", "t", 5, "Threshold of probability")
	identifyCmd.Flags().IntVarP(&o.optParalell, "paralell", "p", 4, "Number of parallel processing")
	identifyCmd.Flags().BoolVarP(&o.optProgress, "pbar", "b", false, "Show progress bar")
	identifyCmd.Flags().IntVarP(&o.optKmer, "kmer", "k", 17, "Length of k-mer")
	identifyCmd.Flags().IntVarP(&o.optSketch, "sketch", "s", 1024, "Sketch size")
}

var identifyCmd = &cobra.Command{
	Use:   "identify",
	Short: "Identifier of plasmid",
	Long:  "Identifier of plasmid using database",
	Run: func(cmd *cobra.Command, args []string) {

		if o.optIn == "" || o.optDB == "" || o.optIdentifyOut == "" {
			fmt.Println("--in, --out and --db are required")
			cmd.Help()
			os.Exit(0)
		}

		inFile := o.optIn
		db := o.optDB
		outFile := o.optIdentifyOut
		NGoRoutines := o.optParalell
		threshold := float64(o.optThreshold) * 0.01
		progress := o.optProgress

		runtime.GOMAXPROCS(NGoRoutines)

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

		binary, err := ioutil.ReadFile(db)
		if err != nil {
			fmt.Println("File reading error", err)
			return
		}

		k := messagePackDecoding(&binary).Kmer
		sketchSize := messagePackDecoding(&binary).SketchSize
		plasmidsMap := messagePackDecoding(&binary).Plasmid

		var bar *pb.ProgressBar
		if progress == true {
			bar = makeBar(&inFile)
			bar.Start()
		}

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

			kmerList := getKmerList(&s, &k, ambiguousDnaComplement)
			minHashValues := make([]*uint64, sketchSize)

			wg := sync.WaitGroup{}
			for i := uint64(0); i < uint64(sketchSize); i++ {
				wg.Add(1)
				go func(i uint64) {
					minValue := uint64((1 << 64) - 1)
					for _, kmer := range *kmerList {
						xxhv := xxhash.Checksum64S(kmer, i)
						if minValue > xxhv {
							minValue = xxhv
						}
					}
					minHashValues[i] = &minValue
					wg.Done()
				}(i)
			}
			wg.Wait()

			var (
				bestHitKey   string
				bestHitValue float64
			)
			similarityMap := make(map[string]float64)

			for refKey, refValue := range plasmidsMap {
				similarity := calcSimilarity(minHashValues, &refValue, 1024)
				similarityMap[refKey] = similarity

				if similarity > bestHitValue {
					bestHitKey = refKey
					bestHitValue = similarity
				}
			}
			fmt.Println(s.Name(), s.Description(), bestHitKey, bestHitValue)

			if bestHitValue >= threshold {
				plasmidSeq.ID = s.Name()
				plasmidSeq.Desc = fmt.Sprintf("Similar to %s (%f)", bestHitKey, bestHitValue)
				plasmidSeq.Seq = s.Slice().(alphabet.Letters)
				_, err := w.Write(plasmidSeq)
				if err != nil {
					log.Printf("Failed to write: %s", err)
				}
			}

			if progress == true {
				bar.Increment()
			}
		}
		if progress == true {
			bar.FinishPrint("Done!")
		}
	},
}

func calcSimilarity(array1 []*uint64, array2 *[]uint64, size float64) float64 {
	var match float64
	for i := 0; i < int(size); i++ {
		if *array1[i] == (*array2)[i] {
			match++
		}
	}
	var bitNum float64 = 64
	similarity := (match / size) - math.Pow(2, -bitNum)/(1-math.Pow(2, -bitNum))

	return similarity
}
