# BEAR - BaseEdit Analyzer
English version | [中文](/README_zh.md)
## Purpose
This program is created to deal with multiple base-edit data in a short time.

One of my colleague once asked me if I can help her to analyse thousands of 
sanger sequence data and show her the edit efficiency of each sample. I was shocked
by the work. She uploaded the ab1 file to a web one by one, which is a tough job.

I know there is a R program, but I know nothing about R and R is very unfriendly
to those known nothing about programing. And R program cannot be packed to
an exe or other executable file. It's very annoying.

So I wrote this program with Python and Qt5. Every user can use it 
even you know nothing about programing! All you need to do is click you mouse
and press "Ctr + C" and "Ctr + V". And the results will be added to **a single EXCEL file**
,which allows you to analyse your data in MicroSoft Excel or WPS.

## Usage
Download the exe file [here](https://github.com/Masterchiefm/BEAR/releases/latest) and execute it.

### Input sample information
Fill the first colum sample name  and the second colum guide is enough. 

If you want the program summarize more information to you, the 3nd, 4th, 5th
colum should be filled. Then it will return the edit efficiency of the unique
position you want.

![](/screenshot.png)

Or you can export the current table, and edit the table in Excel. After edited,
Import the excel table to this program.

### Upload Sanger Sequencing File
For the sequencing file, just select multi samples and drag them to reorganization 
area. Or you can click the 2nd tab to auto fill the table.

After that, click Start! choose a directory and all will be done in minutes.

### Get Results
The results are show here, you can click the "Open Report" button and see the full
report in a local webpage.


In the program, you can see the edit efficiency of each sample and edit efficiency 
at all position.

Most importantly, this table can be export to an Excel file! you don't need
to do more job in data coping!
![](/screenshot2.png)

For the report webpage, you can see the sequencing data and accuracy of your sample, though
most of my user won't open it.

![](/screenshot3.png)


## For Mac Users
Sorry I have no Mac, so I cannot pack it to a mac version. If you have python in your Mac,
download the code and find someone to help you.

Linux users can handle all the difficulties easily, so I have no suggestion.


## Cite This Work
I did't publish it, and I don't know if it can be published in any magazine.
Citing this work in acknowledgement is fine!

If some day I put it to a magazine, I will update this info.
