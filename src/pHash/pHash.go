package main

import (
	"fmt"
	"io"
	"os"

	"./cmd"
)

func main() {
	if err := cmd.RootCmd.Execute(); err != nil {
		if err != io.EOF {
			fmt.Println(err)
			os.Exit(1)
		}
		return
	}
}
