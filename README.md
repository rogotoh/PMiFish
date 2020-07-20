PMiFish
====

## Description
PMiFish: A pipline for MiFish data analysis

## Requirement
* perl 5
* [USEARCH v11](https://www.drive5.com/usearch/): Ultra-fast sequence analysis tool
* [MEGA X](https://www.megasoftware.net/): Molecular Evolutionary Genetics Analysis across computing platforms
* [Gzip](http://gnuwin32.sourceforge.net/packages/gzip.htm) (Windows users require this program. After installing, user need to add manually "GnuWin32\bin" to environment variables)

## Installation
Clone this repository into your local machine  
`git clone https://github.com/rogotoh/PMiFish.git`  

## Usage
1. Put the USEARCH program in Tools folder
2. Put fastq or fastq.gz files generated by Illumina Miseq in Run folder
3. Set each parameters in Setting.txt
4. `perl pmifish.pl`

## Optional analysis
After outputting result files, user can construct phylogenetic trees for each family with database sequences using PA_with_DB.pl
1. Put a setting file (.mao) for MEGA X in Tools folder
2. `perl PA_with_DB.pl`
3. Results files are outputted in 5_2_Summary_Table folder

## Version
2.4.1

## License
[MIT] https://github.com/rogotoh/PMiFish/blob/master/LICENSE
