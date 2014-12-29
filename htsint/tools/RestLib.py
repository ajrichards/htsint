#!/usr/bin/env python
"""
Library of functions to help output data in reStructuredText format
"""

## functions 
def print_rest_table_contents(columns,items,withTrailing=True):
    """
    function needs to be turned into a class
    """

    row = "+"
    head = "+"
    for col in columns:
        row += "-"*col+"-+"
        head += "="*col+"=+"

    if len(columns) != len(items):
        raise Exception("Dimension mismatch %s %s"%(len(columns),len(items)))

    toPrint = "| "

    for i,item in enumerate(items):
        item = str(item)

        if len(item) >= columns[i]:
            raise Exception("col %s not large enough min = %s"%(i,len(item)+2))

        toPrint += item+" "*(columns[i]-len(item))+"| "

    print(toPrint[:-1])

    if withTrailing:
        print(row)
