package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

var (
	o = &Options{}

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
