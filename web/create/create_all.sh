#!/bin/sh
# This gets the environment
. /etc/bashrc


#################################################################
# M. Wobisch 02/16/06
#
# this script creates the whole set of webpages for fastNLO
#
# the basic structure has three levels:
#      homepage -> categories -> subcategories             
#
#  categories are (this numbering scheme will be used later in the code)
#    - 1. news
#    - 2. documentation
#    - 3. calculations
#    - 4. links
#
#  subcategories are so far planned only for: 
#   
#  documentation
#    - concepts: publications, talks
#    - installation: LHAPDF, fastNLO (maybe LHAPDF not needed separately)
#    - user manual: maybe later distinguish user manual and reference manual 
#                          (the latter for the documentation of the table
#                           and the interpolation function, etc.)
#
#  calculations
#    - subcategories for all colliders, if needed Run I, Run II, ...
#    (still need to think how we can list everything within a
#     single subcategory)
#
#################################################################

# -- check version (need version >= 2.0)
echo $BASH_VERSION

# ------- the directories
dir[0]='../'
dir[1]='../news/'
dir[2]='../docs/'
dir[3]='../calculations/'
dir[4]='../form/'
dir[5]='../links/'
readonly dir

#index[0]=



# -------- create index pages for all categories
for i in 0 1 2 3 4 5

do
  cat src_definitions.txt   > tmp-fn.sh
  cat src_header.txt       >> tmp-fn.sh
  cat src_indexbar$i.txt   >> tmp-fn.sh
  cat src_gutter.txt       >> tmp-fn.sh
  cat src_body$i.txt       >> tmp-fn.sh
  cat src_bottom.txt       >> tmp-fn.sh
  echo '_EOF_'             >> tmp-fn.sh

  chmod 744 tmp-fn.sh
  ./tmp-fn.sh               > ${dir[$i]}index.html  
  rm tmp-fn.sh
done

# -------- create empty index page for "create" directory
#          to avoid that people see the source files
echo " empty page ...." > index.html
echo " empty page ...." > ../form/cgi-bin/index.html
echo " empty page ...." > ../form/code/index.html
echo " empty page ...." > ../tables/index.html



# ============== create CGI script =================
cat src_definitions.txt    > tmp-fn1.sh
cat src_cgi02.txt         >> tmp-fn1.sh
chmod 744 tmp-fn1.sh
./tmp-fn1.sh               > src_cgi02html.txt
rm tmp-fn1.sh 

cat src_definitions.txt    > tmp-fn1.sh
cat src_cgi04.txt         >> tmp-fn1.sh
chmod 744 tmp-fn1.sh
./tmp-fn1.sh               > src_cgi04html.txt
rm tmp-fn1.sh 

# ---------
cat src_cgi01.txt          > fn.cgi
cat src_cgi02html.txt     >> fn.cgi
cat src_cgi03.txt         >> fn.cgi
cat src_cgi04html.txt     >> fn.cgi
echo "EOF"                >> fn.cgi

rm src_cgi02html.txt src_cgi04html.txt

# --- copy   cgi-script to target area
chmod 755 fn.cgi
cp fn.cgi ../form/cgi-bin/.
