package cmd

import (
	"fmt"
	"io/ioutil"
	"net/http"
	"os"
	"path"

	"github.com/spf13/cobra"
)

func init() {
	RootCmd.AddCommand(initCmd)
}

var initCmd = &cobra.Command{
	Use:   "init",
	Short: "Download reference plasmid database",
	Long:  "Download reference plasmid database",
	Run: func(cmd *cobra.Command, args []string) {

		url := "https://zenodo.org/record/1991549/files/plasmidDB11062018.phash"

		fmt.Println("Plasmid database is being downloaded...")
		response, err := http.Get(url)
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}

		body, err := ioutil.ReadAll(response.Body)

		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}

		_, filename := path.Split(url)

		file, err := os.OpenFile(filename, os.O_CREATE|os.O_WRONLY, 0666)

		if err != nil {
			fmt.Println(err)
		}

		defer func() {
			file.Close()
		}()

		file.Write(body)
	},
}
