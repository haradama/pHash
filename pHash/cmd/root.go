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
	optIn       string
	optOut      string
	optDB       string
	optParalell int
	optKmer     int
	optSketch   int
}

var (
	mh codec.MsgpackHandle
	o  = &Options{}
)

var RootCmd = &cobra.Command{
	Use:   "pHash",
	Short: "This tool is pretty cool.",
	Long:  "This tool is a great convenience.",
	Run: func(cmd *cobra.Command, args []string) {

	},
}

func init() {
	cobra.OnInitialize()
	RootCmd.AddCommand(versionCmd)
}

var versionCmd = &cobra.Command{
	Use:   "version",
	Short: "Print the version number of go-keisan",
	Long:  `All software has versions. This is go-keisan's`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("pHash v0.1")
	},
}
