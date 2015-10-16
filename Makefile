# ========================================================================
# fl_contour    2.1
# fl_gl_contour 1.9
# geometry      0.7
#
# Copyrigth 2002 by Edmanuel Torres A. (eetorres@gmail.com)
#
# This widget is  free  software;  you  can  redistribute  it and/or
# modify it  under  the  terms  of  the   GNU Library General Public
# License  as  published  by  the  Free  Software Foundation; either
# version 2 of the License,  or  (at your option)  any later version
# or much better FLTK license, wich allow you static linking.
#
# This  library  is  distributed  in the hope that it will be useful,
# but  WITHOUT ANY WARRANTY;  without  even  the  implied warranty of
# MERCHANTABILITY  or FITNESS FOR A PARTICULAR PURPOSE.   See the GNU
# Library General Public License for more details.
#
# You should have  received a copy  of the GNU Library General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA.
#
# Sent me suggestions, modifications or bugs. Don't hesitate to contact
# me for any question, I will be very grateful with your feedbacks.
# Get the last version at http:#eetorres.googlepages.com/fltk_en
#
#========================================================================

include makeinclude

all:   makeinclude
	for dir in $(DIRS); do\
	    echo "==== Compiling in $$dir ====";\
	    (cd $$dir; make || break);\
	done
	
clean: makeinclude
	for dir in $(DIRS); do\
	    echo "==== Cleaning in $$dir ====";\
	    (cd $$dir; make clean || break);\
	done

