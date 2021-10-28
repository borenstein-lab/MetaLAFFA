#!/usr/bin/env python

import argparse
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
from file_handling_lib import *

# Parse command line arguments
parser = argparse.ArgumentParser(
    description="Parses legacy UniRef XML file to create FASTA file and mapping file between UniRef IDs and UniProt IDs.")
parser.add_argument("uniref_xml", help="Legacy UniRef XML file")
parser.add_argument("fasta", help="Output FASTA file location")
parser.add_argument("mapping", help="Output mapping file location")
args = parser.parse_args()

uniref_iterator = ET.iterparse(custom_read(args.uniref_xml))
with open(args.fasta, "w") as fasta, open(args.mapping, "w") as mapping:
    for event, elem in uniref_iterator:
        _, _, elem.tag = elem.tag.rpartition("}")

        # Once we've finished populating the element tree for an entry, write the associated information to the output files
        if elem.tag == "entry":

            entry_id = elem.attrib["id"]
            representative_member = elem.find("representativeMember")

            # FASTA output requires entry ID and protein sequence
            fasta.write("".join([">", entry_id, "\n"]))
            fasta.write("".join([representative_member.find("sequence").text, "\n"]))

            # Mapping output requires entry ID and associated UniProt KB accession IDs, including representative member
            representative_member_id = None
            db_reference = representative_member.find("dbReference")
            if db_reference is not None:
                for elem_property in db_reference.iterfind("property"):
                    if "type" in elem_property.attrib and "value" in elem_property.attrib and elem_property.attrib["type"] == "UniProtKB accession":
                        representative_member_id = elem_property.attrib["value"]
            if representative_member_id is not None:
                mapping.write("".join([entry_id, "\t", representative_member_id, "\n"]))

            for member in elem.iterfind("member"):
                member_db_reference = member.find("dbReference")
                if member_db_reference is not None and "type" in member_db_reference.attrib and member_db_reference.attrib["type"] == "UniProtKB ID":
                    for elem_property in member_db_reference.iterfind("property"):
                        if "type" in elem_property.attrib and "value" in elem_property.attrib and elem_property.attrib["type"] == "UniProtKB accession":
                            mapping.write("".join([entry_id, "\t", elem_property.attrib["value"], "\n"]))
            elem.clear()
