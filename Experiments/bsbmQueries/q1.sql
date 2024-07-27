select ?v0 ?v1 ?v2 ?N1 ?N2 where
{
	?upt <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/ProductType>.
	?v0 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?upt.
	?v0 <http://purl.org/dc/elements/1.1/publisher> ?v1 .
	?v0 <http://purl.org/dc/elements/1.1/date> ?v2 .
	?v0 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/productPropertyNumeric1> ?N1.
	?v0 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/productPropertyNumeric5> ?N2.
}order by (?N1 - ?N2) limit 5
