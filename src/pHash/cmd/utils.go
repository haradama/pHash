package cmd

import (
	"log"
	"unicode/utf8"

	"github.com/ugorji/go/codec"
)

func messagePackEncoding(plasmids *Plasmids) []byte {
	buf := make([]byte, 0, 64)
	err := codec.NewEncoderBytes(&buf, &mh).Encode(*plasmids)
	if err != nil {
		log.Printf("error encoding %v to MessagePack: %v", *plasmids, err)
	}
	return buf
}

func messagePackDecoding(buf *[]byte) Plasmids {
	var plasmids Plasmids
	err := codec.NewDecoderBytes(*buf, &mh).Decode(&plasmids)
	if err != nil {
		log.Printf("error decoding %v to MessagePack: %v", buf, err)
	}
	return plasmids
}

func rev(seq *string) string {
	runes := []rune(*seq)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}

	return string(runes)
}

func getCanonicalKmer(seq string) int {
	num := 0
	for i := 0; i < utf8.RuneCountInString(seq); i++ {
		num += int(seq[i])
	}
	return num
}
