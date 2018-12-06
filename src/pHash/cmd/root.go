package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
	"github.com/ugorji/go/codec"
)

type Plasmids struct {
	Plasmid map[string][]uint64
}

type Options struct {
	optIn          string
	optBuildOut    string
	optIdentifyOut string
	optDB          string
	optParalell    int
	optKmer        int
	optSketch      int
	optThreshold   int
}

var (
	mh codec.MsgpackHandle
	o  = &Options{}
)

var RootCmd = &cobra.Command{
	Use:   "pHash",
	Short: "Software to identify knwon plasmid",
	Long:  "Software to identify knwon plasmid sequence data from metagenome using Minhash",
	Run: func(cmd *cobra.Command, args []string) {

	},
}

func init() {
	cobra.OnInitialize()
	RootCmd.AddCommand(versionCmd)
}

var versionCmd = &cobra.Command{
	Use:   "version",
	Short: "Print the version number of pHash",
	Long:  "Print the version number of pHash",
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("pHash v0.1")
	},
}
