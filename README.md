<img src="./pHash_logo.svg" width="40%" alt="pHash">

pHash is a software to identify known plasmid from metagenomic assembly using MinHash.

## Installation
pHash is available in release page:(https://github.com/haradama/pHash/releases)

## Usage

Please download the plasmid database file on Zenodo: (http://doi.org/10.5281/zenodo.1991549)

```
pHash identify -d plasmidDB11062018.phash -i YOUR_METAGEMOME_FILE
```

## Test
```
sh ./tests/install_test_data.sh
pHash identify -d plasmidDB11062018.phash -i testData.fna
```

## License

[GNU General Public License v3.0](https://github.com/haradama/pHash/blob/master/LICENSE)
