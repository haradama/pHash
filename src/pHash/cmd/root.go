package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
	"github.com/ugorji/go/codec"
)

type (
	Plasmids struct {
		Plasmid map[string][]uint64
	}

	Options struct {
		optIn          string
		optBuildOut    string
		optIdentifyOut string
		optDB          string
		optProgress    bool
		optParalell    int
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

func rev(seq *string) string {
	runes := []rune(*seq)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}

	return string(runes)
}
