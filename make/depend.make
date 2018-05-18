#!/usr/bin/env csh
#
#
#  #################################################################
#  ##                                                             ##
#  ##  depend.make  --  create Makefile dependencies from source  ##
#  ##                   (Generic Unix Version)                    ##
#  ##                                                             ##
#  #################################################################
#
#
foreach x (*.f)
	document 5 $x | grep 'o:'
end
