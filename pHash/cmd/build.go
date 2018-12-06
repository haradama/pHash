package cmd

import (
	"fmt"
	"io"
	"log"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/spf13/cobra"
	"github.com/ugorji/go/codec"

	"./core"
)

func init() {
	RootCmd.AddCommand(buildCmd)
	buildCmd.Flags().StringVarP(&o.optIn, "in", "i", "default", "Input FASTA file")
	buildCmd.Flags().StringVarP(&o.optBuildOut, "out", "o", "reference.phash", "Database")
	buildCmd.Flags().IntVarP(&o.optParalell, "paralell", "p", 4, "Number of parallel processing")
	buildCmd.Flags().IntVarP(&o.optKmer, "kmer", "k", 17, "Length of k-mer")
	buildCmd.Flags().IntVarP(&o.optSketch, "sketch", "s", 1024, "Sketch size")
}

var buildCmd = &cobra.Command{
	Use:   "build",
	Short: "Builder of plasmid database",
	Long:  "Builder of plasmid database using MinHash",
	Run: func(cmd *cobra.Command, args []string) {
		inFile := o.optIn
		outFile := o.optBuildOut
		paralell := o.optParalell
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

		plasmidsMap := map[string][]uint64{}

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
			plasmidsMap[s.Name()] = *minHashValue
		}

		plasmids := Plasmids{
			Plasmid: plasmidsMap,
		}

		file, err := os.Create(outFile)
		if err != nil {
			if err != io.EOF {
				fmt.Println(err)
				os.Exit(1)
			}
		}
		defer file.Close()

		file.Write(messagePackEncoding(plasmids))

	},
}

func messagePackEncoding(plasmids Plasmids) []byte {
	buf := make([]byte, 0, 64)
	err := codec.NewEncoderBytes(&buf, &mh).Encode(plasmids)
	if err != nil {
		log.Printf("error encoding %v to MessagePack: %v", plasmids, err)
	}
	return buf
}
