#!/usr/bin/python

"""Catalytic site search script.

   Automated access to catalytic search server.
   Produces results in comma-separated-value format (standard out).

   Uses curl (via shell process).

   Should work with Python 2.6/7 and 3+.
   """

from __future__ import print_function

import argparse
import HTMLParser
import os
import re
import shlex
import subprocess
import sys
import time

h = HTMLParser.HTMLParser ()
site_url = 'http://catsid.llnl.gov'

# ------------------------------------------------------------------------------
class Debug(object):
    def __init__(self, f):
        self.f = 10 * [False]
        if f:
            for i in f:
                self.f[i] = True


# ------------------------------------------------------------------------------
def main (user_name, password, email, job_id, task_id):

    global do_header_f

    do_header_f = True

    if password:
        if len (password) < 6:
            print ('Password must have at least six characters', file=sys.stderr)
            sys.exit (1)

    if pdb_ids:

        # Get pre-calculated results for these PDBs.
        for pdb_id in pdb_ids:

            # Start search.  -s: silent, except -S: show errors.
            cmd = 'curl -S -s ' \
                  + ' --data-urlencode user_name=%s '          % user_name \
                  + ' --data-urlencode password=%s '           % password \
                  + ' --data-urlencode encrypted_password= ' \
                  + ' --data-urlencode email=%s ' % email \
                  + ' --data-urlencode job_id=%s ' % job_id \
                  + ' --data-urlencode input_protein_name=%s ' % pdb_id \
                  + ' --data-urlencode download_csv_f= ' \
                  + ' %s/find_target' % site_url
            if debug.f[0]:
                print (cmd, file=sys.stderr)

            args = shlex.split (cmd)
            response = subprocess.Popen (args, stderr=sys.stderr, stdout=subprocess.PIPE).communicate ()[0]
            if debug.f[0]:
                print ('pdb-csa matches response:\n' + response, file=sys.stderr)

            print ('Response from server:\n' + response, file=sys.stderr)

        # Done.
        sys.exit (0)
            

    if debug.f[1]:

        # Use previous results.
        if user_name and task_id:
            task_id = user_name + task_id
            ok_f = True
        else :
            if not user_name:
                print ('Need --user_name', file=sys.stderr)

            if not task_id:
                print ('Need --task_id', file=sys.stderr)

            sys.exit (1)

    else:
        ok_f = False
        if not protein_vs_binding_sites:

            # User templates/binding sites vs PDB.  Check that binding site 
            # definitions available.
            if not input_files[0]:
                print ('File of binding-site definitions must be given on command line', file=sys.stderr)
                sys.exit (1)

            if not os.access (define_binding_site, os.R_OK):
                print ('Could not read binding-site definitions file:', define_binding_site, file=sys.stderr)
                sys.exit (1)

            # See if coordinate files given, available.
            for coordinate_file in input_files[1:]:
                if not os.access (coordinate_file, os.R_OK):
                    print ('Could not read binding-site coordinate file:', coordinate_file, file=sys.stderr)
                    sys.exit (1)

            # Upload coordinate files.
            for coordinate_file in input_files[1:]:
                cmd = 'curl -S -s ' \
                      + ' --form user_name=%s ' % user_name \
                      + ' --form password=%s ' % password \
                      + ' --form encrypted_password= ' \
                      + ' --form show_files_flag=0 ' \
                      + ' --form target_search_flag=0 ' \
                      + ' --form bindingsite_protein=bindingsite ' \
                      + ' --form template_file=@%s ' % coordinate_file \
                      + ' %s/user_pdb_file' % site_url
                args = shlex.split (cmd)
                response = subprocess.Popen (args, stderr=sys.stderr, stdout=subprocess.PIPE).communicate ()[0]
                if debug.f[0]:
                    print (cmd, file=sys.stderr)
                    print (response, file=sys.stderr)

                # If no user name, on first pass get assigned user name.
                if not user_name:
                    user_name = user_name_id_task_id (response)

            # Upload binding site file as if entered in textarea -- begins search.
            cmd = 'curl -S -s ' \
                  + ' --data-urlencode user_name=%s ' % user_name \
                  + ' --data-urlencode password=%s ' % password \
                  + ' --data-urlencode encrypted_password= ' \
                  + ' --data-urlencode email=%s ' % email \
                  + ' --data-urlencode job_id=%s ' % job_id \
                  + ' --data-urlencode binding_site_definition@%s ' % define_binding_site \
                  + ' %s/user_binding_sites' % site_url
            args = shlex.split (cmd)
            response = subprocess.Popen (args, stderr=sys.stderr, stdout=subprocess.PIPE).communicate ()[0]
            if debug.f[0]:
                print (cmd, file=sys.stderr)
                print ('PDB search response\n', response, file=sys.stderr)

        else:

            # User proteins vs CSA.
            if not input_files[0]:
                print ('Note: no PDB protein coordinate file(s) given.\nMatching binding sites against previously-uploaded files', file=sys.stderr)
            else:

                # See if target files given, available.
                files_ok_f = True
                for target_file in input_files:
                    if len (target_file) > 50:
                        print ('Protein file name too long (max 50 characters):', target_file, file=sys.stderr)
                        files_ok_f = False

                    if not os.access (target_file, os.R_OK):
                        print ('Could not read protein file:', target_file, file=sys.stderr)
                        files_ok_f = False

                if not files_ok_f:
                    sys.exit (1)

                # Upload files.
                for target_file in input_files:
                    cmd = 'curl -S -s ' \
                          + ' --form user_name=%s ' % user_name \
                          + ' --form password=%s ' % password \
                          + ' --form encrypted_password= ' \
                          + ' --form show_files_flag=0 ' \
                          + ' --form target_search_flag=0 ' \
                          + ' --form bindingsite_protein=protein ' \
                          + ' --form target_file=@%s ' % target_file \
                          + ' %s/user_pdb_file' % site_url
                    args = shlex.split (cmd)
                    response = subprocess.Popen (args, stderr=sys.stderr, stdout=subprocess.PIPE).communicate ()[0]
                    if debug.f[0]:
                        print (cmd, file=sys.stderr)
                        print (response, file=sys.stderr)

                    # If no user name, on first pass get assigned user name.
                    if not user_name:
                        user_name = user_name_id_task_id (response)

            # Start search.
            cmd = 'curl -S -s ' \
                  + ' --data-urlencode user_name=%s ' % user_name \
                  + ' --data-urlencode password=%s ' % password \
                  + ' --data-urlencode encrypted_password= ' \
                  + ' --data-urlencode email=%s ' % email \
                  + ' --data-urlencode job_id=%s ' % job_id \
                  + ' --data-urlencode show_files_flag=0 ' \
                  + ' --data-urlencode target_search_flag=1 ' \
                  + ' --data-urlencode bindingsite_protein=protein ' \
                  + ' %s/user_pdb_file' % site_url
            args = shlex.split (cmd)
            response = subprocess.Popen (args, stderr=sys.stderr, stdout=subprocess.PIPE).communicate ()[0]
            if debug.f[0]:
                print (cmd, file=sys.stderr)
                print ('protein search response:\n', response, file=sys.stderr)

        # Retrieve user name and task ID from hidden form fields.
        #task_id = user_name_id_task_id (response)

        # Wait for result.
        #time.sleep (10)

    # Response should be message on results location.
    print (response)


    # If result not available, wait again.
    #for i in ['20']:

    #   # Turn globbing off - assigned user names (which begin with '}') are
    #   # a problem.
    #   cmd = 'curl -S -s --globoff "%s/search_results/%s/"' \
    #         % (site_url, task_id)
    #   if debug.f[0]:
    #       print (cmd, file=sys.stderr)

    #   args = shlex.split (cmd)
    #   response = subprocess.Popen (args, stderr=sys.stderr, stdout=subprocess.PIPE).communicate ()[0]
    #   match = re.search ('results_table', response)
    #   if match:
    #       ok_f = True
    #       break

    #   time.sleep (20)

    #if not ok_f:
    #   print ('Did not get search results. Here is the last response from the server.\n', file=sys.stderr)
    #   print (response, file=sys.stderr)
    #else:
    #   if debug.f[0]:
    #       print ('search results\n' + response, file=sys.stderr)

    #   # Parse results from table to csv.
    #   parse_table_to_csv (response, match)


# ------------------------------------------------------------------------------
def parse_table_to_csv (response, match, pdb_id=''):
    global do_header_f

    start = match.start (0)

    # Make each table row a list element.
    rows = re.findall (r'<tr[^>]*>(.*?)</tr>', response[start:], re.DOTALL)
    if debug.f[0]:
        print ('[parse_table_to_csv] len (rows)', len (rows), file=sys.stderr)

    for i_row, row in enumerate (rows):

        # Don't do header if already done.
        if i_row == 0 and not do_header_f:
            continue

        # Stop if we've come to the footer (which is also a table).
        if 'Lawrence Livermore National Laboratory' in row:
            break

        # Make each table item a list element.
        items = re.findall (r'<t[dh][^>]*>(.*?)</t[dh]>', row, re.DOTALL)

        # Skip if only two items (additional rows just for alternate UniProt
        # EC numbers).
        if len (items) == 2:
            continue

        new_items = []
        if pdb_id:
            if i_row == 0:
                if do_header_f:
                    new_items.append ('pdb_id')
            else:
                new_items.append (pdb_id)

        # Skip the first, "view" item.
        for item in items[1:]:
            new_item = item.strip ()

            # Delete any tags, newlines.
            new_item = re.sub (r'<[^>]*>', '', new_item)
            new_item = re.sub ('\n', ' ', new_item)

            # Convert multiple spaces to single space.
            new_item = re.sub (' +', ' ', new_item)

            # Convert any HTML entities (like &#39 for single quote or 
            # "prime") back to unicode.
            new_item = h.unescape (new_item)

            new_items.append (new_item.encode ('utf-8'))


        print (','.join (new_items))

        # If have come to a score below cutoff, stop.
        if pdb_id and i_row > 0 and float (items[1]) < cutoff_score:
            break


    # Flag that have written header.
    do_header_f = False


# ------------------------------------------------------------------------------
def user_name_id_task_id (response):
    match = re.search (r'name="user_name" +value="([^"]+)"', response)
    user_name_id = match.group (1)
    user_task_id = user_name_id

    match = re.search (r'name="task_id" +value="(\w+)"', response)
    if match:
        task_id = match.group (1)
        user_task_id += '__' + task_id

    if debug.f[0]:
        print ('user_task_id', user_task_id, file=sys.stderr)

    return user_task_id


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

cutoff_score = 0.001

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Catalytic site identification script -- command-line access to server.  Produces comma-separated values (standard out).  Please be patient - may take a couple of minutes.')

    parser.add_argument('--email', metavar='ADDRESS', help='Email address (required)', required='true')
    parser.add_argument('--job_id', metavar='JOB_ID', help='Name for this batch job (required)', required='true')

    parser.add_argument('--protein_vs_binding_sites', help='Search for binding sites in user-supplied proteins (default: search PDB for matches to user-supplied binding sites).', action='store_true', default=False)
    parser.add_argument('--pdb_ids', nargs='+', metavar='PDB_ID', help='List of PDB entries for which matches to catalytic sites will be shown.')
    parser.add_argument('--cutoff_score', metavar='CUTOFF', help='Score below which results will not be printed (default: %s )' % cutoff_score)
    parser.add_argument('--user_name', metavar='USER_NAME', help='Optional name under which uploaded files and results will be saved on server.')
    parser.add_argument('--password', metavar='PASSWORD', help='Optional password for results (default: no password)')
    parser.add_argument('--task_id', metavar='TASK_ID', help='For debug 1, from results url: search_results/<user_name><task_id>/ (for example, joeuser__67)')
    parser.add_argument('--debug', type=int, nargs='+', metavar='DEBUG_NUM',  help='Debug numbers to turn on (0: response from uploads; 1: use previous results -- task_id required)')

    parser.add_argument('define_binding_site', nargs='?', metavar='DEFINE_BINDING_SITE', help='File of binding site residue definitions (required for default search of PDB with user-supplied binding sites)')
    parser.add_argument('coordinate_files', nargs='*', metavar='COORDINATE_FILE', help='File(s) of binding site residue coordinates (default) or protein coordinates (if --protein_vs_binding_sites)')

    args = parser.parse_args()

    email                    = args.email
    if not re.match("^[a-zA-Z0-9._%-]+@[a-zA-Z0-9._%-]+\.[a-zA-Z]{2,6}$", email):
        print ('Email address not valid: ', email)
        sys.exit (1)

    job_id                   = args.job_id
    protein_vs_binding_sites = args.protein_vs_binding_sites
    pdb_ids                  = args.pdb_ids
    if protein_vs_binding_sites and pdb_ids:
        print ('Cannot do both --protein_vs_binding_sites and --pdb_ids', file=sys.stderr)
        sys.exit (1)

    user_name                = args.user_name
    if not user_name:
        user_name = ''

    password                 = args.password
    if not password:
        password = ''

    cutoff_score             = args.cutoff_score
    task_id                  = args.task_id

    define_binding_site      = args.define_binding_site
    coordinate_files         = args.coordinate_files

    input_files = [define_binding_site] + coordinate_files

    debug       = Debug(args.debug)

    main (user_name, password, email, job_id, task_id)

