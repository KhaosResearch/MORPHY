#!/bin/bash 
red='\033[0;31m'
blue='\033[1;34m'
green='\033[0;32m'
NC='\033[0m' 
bold='\033[1m'
cyan='\033[0;36m'


SCRIPTPATH=$PWD

echo -e "${cyan}----------- ${bold}MORPHY Installer${cyan} ------------${NC}"
echo -e "${green}${bold}1.- Install Bio++ (Bpp-Core, Bpp-Seq and Bpp-Phyl Libraries)${NC}"

bpp_dir=$SCRIPTPATH/lib/Bpp
echo -e "${blue}${bold}Change to Bpp Directory: $bpp_dir${NC}"
cd $bpp_dir

for d in bpp-core-2.1.0 bpp-seq-2.1.0 bpp-phyl-2.1.0; do
  if [[ ! -e $d ]]; then
    echo -e "${red}${bold} Skip $d (does not exist).{NC}"
    continue
  fi
  
  echo -e "${blue}${bold}Installing  $d...${NC}"
  cd $d/
  cmake -D CMAKE_INSTALL_PREFIX=/usr/local -D BUILD_TESTING=FALSE ./
 
  if make install; then
	 echo -e "${cyan}${bold}$d Installed.!!!${NC}"
  else
	echo -e "${red}${bold}Compilation error in '$d/'. Abort${NC}"       
        cd ../..
        break
  fi
  cd ..
done

echo -e "${green}${bold}2.- Install Phylogenetic Likelihood Library (PLL)${NC}"

cd $SCRIPTPATH/lib/Pll

echo -e "${blue}${bold}autoreconf${NC}"
autoreconf -fvi 

echo -e "${blue}${bold}./configure${NC}"
./configure 

echo -e "${blue}${bold}Installing Pll${NC}"
if make install; then
	 echo -e "${cyan}${bold}Pll Installed.!!!${NC}"
  else
	echo -e "${red}${bold}Compilation error Pll.${NC}"       
  fi

echo -e "${green}${bold}3.- Install MORPHY${NC}"
cd $SCRIPTPATH/

echo -e "${blue}${bold}Installing MORPHY${NC}"
if make; then
	echo -e "${cyan}${bold}MORPHY is Installed${NC}"
  else
	echo -e "${red}${bold}Compilation error MORPHY.${NC}"       
fi







