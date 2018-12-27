package cmd

import (
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
	pb "gopkg.in/cheggaaa/pb.v1"
)

func init() {
	RootCmd.AddCommand(buildCmd)
	buildCmd.Flags().StringVarP(&o.optIn, "in", "i", "default", "Input FASTA file")
	buildCmd.Flags().StringVarP(&o.optBuildOut, "out", "o", "reference.phash", "Database")
	buildCmd.Flags().IntVarP(&o.optParalell, "paralell", "p", 4, "Number of parallel processing")
	buildCmd.Flags().BoolVarP(&o.optProgress, "pbar", "b", false, "Show progress bar")
	buildCmd.Flags().IntVarP(&o.optKmer, "kmer", "k", 17, "Length of k-mer")
	buildCmd.Flags().IntVarP(&o.optSketch, "sketch", "s", 1024, "Sketch size")
}

var buildCmd = &cobra.Command{
	Use:   "build",
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
		NGoRoutines := o.optParalell
		k := o.optKmer
		sketchSize := o.optSketch
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

		plasmidsMap := map[string][]uint64{}

		var bar *pb.ProgressBar
		if progress == true {
			bar = makeBar(&inFile)
			bar.Start()
		}

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
			minHashValues := make([]uint64, sketchSize)

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
					minHashValues[i] = minValue
					wg.Done()
				}(i)
			}
			wg.Wait()

			plasmidsMap[s.Name()] = minHashValues

			if progress == true {
				bar.Increment()
			}
		}

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

		if progress == true {
			bar.FinishPrint("Done!")
		}
	},
}
