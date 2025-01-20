# pyrodiversity
Code to calculate pyrodiversity according to [Steel et al. 2021](https://royalsocietypublishing.org/eprint/RCRFTQ4V7Z2C6SFWWXG3/full). 

Includes demos for [prepping fire perimeters](https://github.com/zacksteel/pyrodiversity/blob/master/code/PerimeterPrep.md) for burn severity estimation in Google Earth Engine, and for [calculating pyrodiversity](https://github.com/zacksteel/pyrodiversity/blob/master/code/YosemiteDemo.md) at the watershed-scale using Yosemite National Park. Pyrodiversity can also be quantified around sample points or plots. A demo for doing this is coming soon.  

If using this code please cite:  
Z. L. Steel, B. M. Collins, D. B. Sapsis, and S. L. Stephens, “Quantifying pyrodiversity and its drivers,” Proceedings of the Royal Society B: Biological Sciences, vol. 288, no. 1948, Apr. 2021.

2022-04-26 Update:  
Surface generation functions and global_fd() (formerly known as pyrodiv_calc() ) have been updated to improve functionality and increase speed. Please let me know if you find any bugs. 

2025-01-19:
Periodic bug fixes, error checking, and improved messaging in surface and *fd functions.
