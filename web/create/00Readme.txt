

******************************************************************
# M. Wobisch 02/16/06
#
# the fastNLO webpages
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
******************************************************************


* never edit one of the final webpages! They are created 
  automatically by the script below.


******************************************************************
**** M Wobisch 02/16/06
****
**** how to create the fastNLO webpages
****
******************************************************************

* go to the directory in which the file index.html for the homepage
  should reside. 
  - create subdirectores "news", "docs", calculations", "links"
  - create another subdirectory "create" and copy all the
    files from this package into  create/.

* edit the file "src_definitions.txt"
  to set the absolute path of your webpages
  (and the desired properties - colors and the size)

* run   ./create_all.sh
  this creates all webpages in the top-level directory as well 
  as in "news", "docs", calculations", "links"

* voila!


******************************************************************

There are six elements from which the homepage and 
each catagory main-page are built:

src_definitions.txt    universal - contains the fundamental 
                       definitions of colors etc.
src_header.txt         universal - defines the header
src_indexbar#.txt      category-specific  (# is from 0-4)
                       defines the indexbar on the left
src_gutter.txt         universal - just a small separator
                       in here the height is defined
src_body#.txt          category-specific  (# is from 0-4)
                       defines the main area ("body" - # is from 0-4)
src_bottom.txt         universal - the bottom area
              
remember: the top-header and the date at the bottom
          are added by CEDAR

******************************************************************

The subcategories .... to be implemented

******************************************************************
