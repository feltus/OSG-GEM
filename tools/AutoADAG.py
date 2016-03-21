# #
#  Copyright 2007-2012 University Of Southern California
#
#  Licensed under the Apache License, Version 2.0 (the 'License');
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing,
#  software distributed under the License is distributed on an 'AS IS' BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# #

__author__ = 'Rajiv Mayani'

import logging

try:
    from Pegasus.DAX3 import ADAG, Job, File, Executable, PFN, Link, When, DuplicateError
except ImportError, e:
    logging.error('Include Pegasus Python libraries in your PYTHONPATH')


class AutoADAG(object, ADAG):
    """
    Automatically determine the dependencies between jobs based on the file usages.
    All jobs consuming a file F depend on the singular job that produces that file.
    """
    def __init__(self, name, count=None, index=None):
        ADAG.__init__(self, name, count, index)

    def writeXML(self, out):

        mapping = {}

        def addOutput(job, file_obj):

            if file_obj:
                file_obj = file_obj.name

                if file_obj not in mapping:
                    mapping[file_obj] = (set(), set())

                mapping[file_obj][1].add(job)

        # Automatically determine dependencies

        # Traverse each job
        for job_id, job in self.jobs.iteritems():
            file_used = job.used

            # If job produces to stdout, identify it as an output file
            addOutput(job, job.stdout)
            # If job produces to stderr, identify it as an output file
            addOutput(job, job.stderr)

            # If job consumes from stdin, identify it as an input file
            if job.stdin:
                if job.stdin.name not in mapping:
                    mapping[job.stdin.name] = (set(), set())

                mapping[job.stdin.name][0].add(job)


            for file in file_used:

                if file.name not in mapping:
                    mapping[file.name] = (set(), set())

                if file.link == Link.INPUT:
                    mapping[file.name][0].add(job)
                else:
                    mapping[file.name][1].add(job)

        for file_name, io in mapping.iteritems():

            # Go through the mapping and for each file add dependencies between the
            # job producing a file and the jobs consuming the file
            inputs = io[0]

            if len(io[1]) > 0:
                output = io[1].pop()

                for input in inputs:
                    try:
                        self.depends(parent=output, child=input)
                    except DuplicateError:
                        pass

        super(AutoADAG, self).writeXML(out)
