#!/usr/bin/python

import re
import sys

g_processed_msg = '\x3c\x21\x2d\x2d\x20\x28\xc2\xaf\x60\xc2\xb7\x2e\x5f\x2e\xc2\xb7\xc2\xb7\xc2\xb8\x2e\x2d\x7e\x2a\xc2\xb4\xc2\xa8\xc2\xaf\xc2\xa8\x60\x2a\xc2\xb7\x7e\x2d\x2e\x2c\x2d\x28\x20\x50\x4f\x53\x54\x2d\x50\x52\x4f\x43\x45\x53\x53\x45\x44\x20\x42\x59\x20\x50\x4c\x4c\x20\x44\x4f\x58\x59\x47\x45\x4e\x20\x50\x52\x4f\x43\x45\x53\x53\x4f\x52\x20\x29\x2d\x2c\x2e\x2d\x7e\x2a\xc2\xb4\xc2\xa8\xc2\xaf\xc2\xa8\x60\x2a\xc2\xb7\x7e\x2d\x2e\xc2\xb8\xc2\xb7\xc2\xb7\x2e\x5f\x2e\xc2\xb7\xc2\xb4\xc2\xaf\x29\x20\x20\x2d\x2d\x3e\n'

g_path_prefix = '/var/www/test/html/'
#  files to process
g_html_groups = ('group__evaluateLikelihoodGroup.html', 'group__rearrangementGroup.html', 'group__alignmentGroup.html', 'group__newickParseGroup.html', 'group__parsePartitionFileGroup.html')
g_html_files  = ('newviewGenericSpecial_8c.html', 'optimizeModel_8c.html', 'genericParallelization_8c.html', 'utils_8c.html', 'evaluateGenericSpecial_8c.html', 'makenewzGenericSpecial_8c.html')

g_memberdecls_table_string = ('table', 'class', 'memberdecls')
g_bStatic = 0
g_bNormal = 1
g_static_funtable_start = '<table class="memberdecls"><tr class="heading"><td colspan="2"><h2 class="groupheader"><a name="func-members"></a>Static functions</h2></td></tr>\n'
g_static_funtable_end   = '</table>\n'

def locate_memberdecls_index (content):
  dfa = re.compile (".*".join (g_memberdecls_table_string))
  for i in range(len(content)):
    if dfa.search (content[i]) is not None: break

  return i
   

def readFile (fname):
  with open(fname) as f:
    content = f.readlines()

  return content

def split_html_content (content):
  i = locate_memberdecls_index (content)

  # nothing found or found in the last line, quit
  if i == len (content) - 1: return (content, [], [])

  dfa = re.compile ("^</table>")
  # get the ending </table> tag
  for j in range (i + 3, len(content)):
    if dfa.search (content[j]) is not None: break
  
  header    = content[:i + 3]
  functions = content[i + 3:j]
  footer    = content[j:]

  return (header, functions, footer)

def split_static_functions (functions):
  dfa_fun   = re.compile ("^<tr[^>]*class[^>]*memitem.*><td((?!td>).)*static", re.DOTALL)
  dfa_desc  = re.compile ("^<tr[^>]*class[^>]*memdesc[^>]*>")
  dfa_trend = re.compile ("</tr>\n$");
  static = []
  normal = []

  where = 0

  i = 0
  while (i < len(functions)):
    #while len (functions[i].strip()) == 0: i = i + 1
    
    # we assume lines end with </tr> -- read function declaration line
    line = ""
    while dfa_trend.search (line) is None:
      line = line + functions[i]
      i = i + 1

    if dfa_fun.search (line) is not None:
      where = g_bStatic
      static.append (line)
    else:
      normal.append (line)
      where = g_bNormal


    # we assume lines end with </tr> -- read optional @brief oneliner
    line = ""
    while dfa_trend.search (line) is None:
      line = line + functions[i]
      i = i + 1

    # check if there is a @brief description and append it to the appropriate list
    if dfa_desc.search (line) is not None:
      if where == g_bStatic:
        static.append (line)
      else:
        normal.append (line)

      #if there was a @brief oneliner, now read the horizontal bar
      line = ""
      while dfa_trend.search (line) is None:
        line = line + functions[i]
        i = i + 1
      #i = i + 1

    # store the horizontal bar
    if where == g_bStatic:
      static.append (line)
    else:
      normal.append (line)

  return static, normal

def overwrite_file (fname, header, static, normal, footer):
  
  f = open (fname, 'w')

  f.write(g_processed_msg)

  normalExist = True


  for line in header:
    f.write(line)


  if len(normal) != 0:
    for line in normal:
      f.write(line)
  else:
    normalExist = False

  if normalExist == True and len(static) != 0:
    f.write(g_static_funtable_start)
    for line in static:
      f.write(line)
    f.write(g_static_funtable_end)
  else:
    if normalExist == False:
      for line in static:
        f.write(line)

  for line in footer:
    f.write(line)
  
  f.close()


def process_static_functions (fname):
  """Separate static from normal functions in file fname and overwrite the file"""
  content = readFile (fname)

  # check if file was already processed 
  if len(content) > 0:
    if content[0] == g_processed_msg:
      print 'File ' + fname + ' already processed...'
      return


  header, functions, footer = split_html_content (content)

  # no functions found? quit
  if not functions: return
  
  static, normal = split_static_functions (functions)

  overwrite_file (fname, header, static, normal, footer)


if __name__ == "__main__":
  html_files = g_html_groups + g_html_files
  for f in html_files:
    process_static_functions (g_path_prefix + f)
