def vmet_logo() {
figlet =
"""\
                           _
__   __     _ __ ___   ___| |_
\\ \\ / /____| '_ ` _ \\ / _ \\ __|
 \\ V /_____| | | | | |  __/ |_
  \\_/      |_| |_| |_|\\___|\\__|
"""

version = "v${workflow.manifest.version}".center(31)

return \
"""\
${figlet}
${version}\
"""
}
