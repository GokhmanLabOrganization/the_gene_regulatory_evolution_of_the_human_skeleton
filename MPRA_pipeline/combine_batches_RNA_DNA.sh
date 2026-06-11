cells=chondrocytes
library=L1
adaptor=a1

for f in /home/labs/davidgo/Collaboration/humanMPRA/$cells/$library$adaptor/input/DNA_RNA/batch1/*; do name=$(basename $f ); cat /home/labs/davidgo/Collaboration/humanMPRA/$cells/$library$adaptor/input/DNA_RNA/batch1/$name /home/labs/davidgo/Collaboration/humanMPRA/$cells/$library$adaptor/input/DNA_RNA/batch2/$name > /home/labs/davidgo/Collaboration/humanMPRA/$cells/$library$adaptor/input/DNA_RNA/$name; echo $name; done