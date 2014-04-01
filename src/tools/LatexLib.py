#!/usr/bin/env python
import os,sys,time,re
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

## create doc beginning
class LatexReportCreator():
    def __init__(self,docName,docType='report',docDir="report",verbose=False,extraPackages=[],generatedBy="Python"):
        self.docName = docName
        self.docType = docType
        self.docDir = docDir
        self.verbose = verbose
        self.generatedBy = generatedBy
        self.extraPackages = extraPackages
        self.fontSizeList = ["tiny","scriptsize","footnotesize","small","normalsize","large","Large","LARGE","huge","Huge"] 

        ## error check
        if type(self.extraPackages) != type([]):
            print "INPUT ERROR: extra packages must be in list form"

        if os.path.isdir(self.docDir) == False:
            if self.verbose == True:
                print "...making report directory"
            os.mkdir(self.docDir)
        
        self.fid = open(os.path.join(self.docDir,self.docName + ".tex"),'w')

    def write_beginning(self,title,author,shortTitle='report',shortAuthor='',toc=False):
        self.fid.write("\documentclass[letterpaper,10pt]{%s}\n"%self.docType)
        self.fid.write("\usepackage{pgf,graphics,color,hyperref,fullpage,natbib,algorithm,algorithmic,amsmath,ulem,bm,morefloats}\n")
        for package in self.extraPackages:
            self.fid.write("\usepackage{%s}\n"%package)
        self.fid.write("\definecolor{darkblue}{rgb}{0.0,0.0,0.50}\n")
        self.fid.write("\definecolor{darkgreen}{rgb}{0.0,0.35,0.0}\n")
        self.fid.write("\hypersetup{colorlinks=true, linkcolor=darkblue, citecolor=darkblue, urlcolor=darkblue}\n")
        self.fid.write("\hypersetup{pdfauthor={%s}, pdftitle={%s}}\n"%(shortAuthor,shortTitle))
        self.fid.write("\usepackage{titlesec}\\titleformat{\chapter}{\centering\\normalfont\Large\\bfseries}{\\thechapter}{1em}{}\n")

        ## beginning materials
        self.fid.write("\\begin{document}\n")
        self.fid.write("%%%%%\n")
        self.fid.write("\\begin{flushleft}\n")
        self.fid.write("\\textbf{%s}\\\ \n"%title)
        self.fid.write("Created by: %s\\\ \n"%author)
        self.fid.write("Created on: %s\n"%time.asctime())
        self.fid.write("\end{flushleft}\n")
        self.fid.write("\medskip\hrule height 1pt\n")
        self.fid.write("\\vspace{7pt}\n")

        if toc == True:
            self.fid.write("\\tableofcontents\n")

    ## convenience functions
    def begin(self,item):
        if item in self.fontSizeList:
            self.fid.write("\%s\n"%item)
        else:
            self.fid.write("\\begin{%s}\n"%item)

    def end(self,item):
        if item in self.fontSizeList:
            self.fid.write("\\normalsize\n")
        else:
            self.fid.write("\end{%s}\n"%item)

    def section(self,sectionName,sectionType="section",numbers=True,toc=False):
        '''
        toc - used to manually add a section to the table of contents -- must also call toc = True in 'write_beginning'

        '''

        if numbers == True:
            self.fid.write("\%s{%s}\n"%(sectionType,sectionName))
        else:
            self.fid.write("\%s*{%s}\n"%(sectionType,sectionName))

        if toc == True:
            self.fid.write("\\addcontentsline{toc}{%s}{%s}\n"%(sectionType,sectionName))

    ## main functions
    def write_paragraph(self,body):
        self.fid.write(body+"\n\n")

    def write_end(self):
        self.fid.write("\n\end{document}")
        self.fid.close()

    def compile(self):
        currentDir = os.getcwd()
        os.chdir(self.docDir)
        os.system("pdflatex %s"%self.docName+".tex")
        os.chdir(currentDir)

    def include_figure(self,figPath,caption=None,label=None,scale=0.65):
        if os.path.isfile(figPath) == False:
            print "ERROR: could not include figure does not exist", figPath

        self.fid.write("\\begin{figure}[!ht]\n")
        self.begin("center")
        self.fid.write("\includegraphics[scale = %s]{%s}\n"%(scale,figPath))
        if caption != None:
            self.fid.write("\caption{%s}\n"%caption)
        if label != None:
            self.fid.write("\label{fig:%s}\n"%label)
        self.end("center")
        self.end("figure")

    def include_table(self,rowDict,colHeader=None,customHeader=None,justifications=None,caption=None,label=None,fontSize='normalsize',ordered=False):
        '''
        rowDict - which has row header elements as keys and lists to fill out the cols as values
        colHeader - are the header elements for the columns
        '''
        
        ## if keys are only numeric
        numericKeys = False
        if not re.search('\D',str(rowDict.keys()[0])): 
            numericKeys = True

        ## error checking
        if type(rowDict) != type({}):
            print "ERROR: rowDict for include_table must be of type dict"
        if type(rowDict[rowDict.keys()[0]]) != type([]):
            print "ERROR: the values in rowDict for included_table must be of type list"

        numCols = len(rowDict[rowDict.keys()[0]])
        if justifications == None and numericKeys == False:
            justifications = '|l|'+'c'*numCols+"|"
        elif justifications == None and numericKeys == True:
            justifications = '|l|'+'c'*(numCols-1)+"|"
            
        self.fid.write('\\begin{table}[!h]\n')
        self.begin('center')
        self.begin(fontSize)
        self.fid.write("\\begin{tabular}{%s}\n"%justifications)
        self.fid.write("\\hline\n")

        if customHeader != None:
            self.fid.write(customHeader)
            self.fid.write("\\hline\n")
        
        if colHeader != None:
            header = "".join([i + "&" for i in colHeader])[:-1]+"\\\ \n"
            self.fid.write(header)
            self.fid.write("\\hline\n")

        if numericKeys == True:      
            sk = [int(k) for k in rowDict.keys()]
            sk.sort()
            for key in sk:
                row = rowDict[str(key)]
                if row[0] == "hline":
                    self.fid.write("\\hline\n")
                else:
                    row = "".join([i + "&" for i in row])[:-1]+"\\\ \n"
                    self.fid.write(row)
        elif type(ordered) == type([]):
            for rowHead in ordered:
                row = rowDict[rowHead]
                if rowHead == 'hline':
                    self.fid.write("\\hline\n")
                else:
                    row = rowHead + "&" + "".join([i + "&" for i in row])[:-1]+"\\\ \n"
                    self.fid.write(row)
        elif ordered == True:
            it = iter(sorted(rowDict.iteritems()))
            for rowHead, row in it:
                if row[0] == 'hline':
                    self.fid.write("\\hline\n")
                else:
                    row = rowHead + "&" + "".join([i + "&" for i in row])[:-1]+"\\\ \n"
                    self.fid.write(row)
        else:
            for rowHead, row in rowDict.iteritems():
                if row[0] == 'hline':
                    self.fid.write("\\hline\n")
                else:
                    row = rowHead + "&" + "".join([i + "&" for i in row])[:-1]+"\\\ \n"
                    self.fid.write(row)

        self.fid.write("\\hline\n")
        self.end('tabular')
        self.end(fontSize)
        self.end('center')
        if caption != None:
            self.fid.write("\caption{%s}\n"%caption)
        if label != None:
            self.fid.write("\label{tab:%s}\n"%label)
        self.end('table')

## run a small test 
if __name__ == '__main__':

    ## create a preambble
    report = LatexReportCreator('test',verbose=True)
    report.write_beginning("Test Report Title","Ender Wiggin",toc=True)
    
    ## create the body
    report.section("About",numbers=False,toc=True)
    
    report.write_paragraph("this is a paragraph blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah "+ 
                           "blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah")
    report.write_paragraph("this is a paragraph blah blah blah blah blah blah blah blah blah blah blah blah "+ 
                           "blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah")
   
    report.section("Results",numbers=False,toc=True)
    report.write_paragraph("this is a paragraph blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah "+
                           "blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah")
    report.write_paragraph("this is a paragraph blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah "+
                           "blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah")
   
    colHeader = ["I", "II", "III", "IV"] 
    rowDict = {"a":["1","2","3"],
               "b":["4","5","6"],
               "c":["7","8","9"]   
               }

    caption = "\\textbf{an example table}. this is a table caption"
    label = "example-table"
    report.include_table(rowDict,colHeader,caption=caption,label=label)

    ## create a test figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n = 5000
    alpha,beta = 2.0,3.0
    samples = np.random.beta(alpha,beta,size=n)
    p1 =ax.hist(samples,normed=1,alpha=0.7)
    pdfX = np.linspace(0,1,300)
    pdfY = stats.beta.pdf(pdfX,alpha,beta)
    p2 = ax.plot(pdfX,pdfY,color='orange',linewidth=2.5)
    ax.set_xlim((0,1))

    ax.set_title(r"Beta distribution: $a$ = %s, $b$ = %s, $n$ = %s"%(alpha,beta,n))
    figPath = os.path.join(report.docDir,"examplefig.pdf")
    ax.legend([p2],["pdf"])
    plt.savefig(figPath)

    ## add in test figure
    caption = "\\textbf{The beta distribution}. this is a figure caption"
    label = "beta-distn"
    report.include_figure("examplefig.pdf",caption=caption,label=label,scale=0.5)

    ## add another paragraph and close report
    report.write_paragraph("This is an example of how to make a reference to Figure~\\ref{fig:beta-distn} blah blah blah blah blah "+
                           "blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah")

    report.write_end()
    report.compile()
