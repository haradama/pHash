package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
	"github.com/ugorji/go/codec"
)

type (
	Plasmids struct {
		SketchSize int
		Kmer       int
		Plasmid    map[string][]uint64
	}

	Options struct {
		optIn          string
		optBuildOut    string
		optIdentifyOut string
		optDB          string
		optKmer        int
		optSketch      int
		optThreshold   int
	}
)

var (
	mh codec.MsgpackHandle
	o  = &Options{}

	RootCmd = &cobra.Command{
		Use:   "pHash",
		Short: "Software to identify knwon plasmid",
		Long:  "Software to identify knwon plasmid from metagenome using Minhash",
		Run: func(cmd *cobra.Command, args []string) {
			cmd.Help()
		},
	}

	versionCmd = &cobra.Command{
		Use:   "version",
		Short: "Print the version number of pHash",
		Long:  "Print the version number of pHash",
		Run: func(cmd *cobra.Command, args []string) {
			fmt.Println("pHash v0.2")
		},
	}
)

func init() {
	cobra.OnInitialize()
	RootCmd.AddCommand(versionCmd)
}
