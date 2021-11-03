#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def cowsay(message) {
    messagelines = message.split('\n')
    nlines = messagelines.length
    linelength = 0
    messagelines.each{ l -> if ( l.length() > linelength ) { linelength = l.length() } }
    paddinglength = linelength + 2

    if ( nlines == 1 ) {
        balloon =
            """ ${"_"*paddinglength}
            < ${message} >
             ${"-"*paddinglength}"""
    }
    else {
        balloon =
            """ ${"_"*paddinglength}
/ ${messagelines[0].padRight(linelength)} \\"""
        for (int i=1;i<(nlines-1);i++) {
            balloon =
                """${balloon}
| ${messagelines[i].padRight(linelength)} |"""
        }
    balloon =
    """${balloon}
\\ ${messagelines[nlines-1].padRight(linelength)} /
 ${"-"*paddinglength}"""
    }

    cow =
    """    \\   ^__^
         \\  (oo)\\_______
            (__)\\       )\\/\\
                ||----w |
                ||     ||"""

    cowput =
        """${balloon}
    ${cow}
"""

    log.info cowput
}
