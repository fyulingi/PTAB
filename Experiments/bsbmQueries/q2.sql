select ?v0 ?v1 ?v2 ?v3 ?N1 ?N2 ?N3  where
{
	?v0 <http://www.w3.org/2000/01/rdf-schema#label> ?v1 .
	?v0 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?v2 .
	?v3 <http://www.w3.org/2000/01/rdf-schema#label> ?v1 .
	?v3 <http://purl.org/dc/elements/1.1/publisher> <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/instances/StandardizationInstitution1> .
	?v0 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/productPropertyNumeric1> ?N1.
	?v0 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/productPropertyNumeric2> ?N2.
	?v0 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/productPropertyNumeric3> ?N3.
}order by (-?N1 - ?N2 -?N3) limit 5
