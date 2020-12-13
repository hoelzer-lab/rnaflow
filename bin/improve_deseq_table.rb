#!/usr/bin/env ruby

######################
## input is a DEGs table in csv format
## and a hash file between ID and real gene name
## build new csv file with extended information

tmp_out = File.open(ARGV[0],'w')
csv = File.open(ARGV[1],'r')
hash = File.open(ARGV[2],'r')

h = {}
hash.each do |line|
  s = line.split("\t")
  id = s[0]; name = s[1]; type = s[2].chomp; chr = s[3].chomp; start = s[4].chomp; stop = s[5].chomp; strand = s[6].chomp; desc = s[7].gsub(',',';').chomp
  h[id] = [name, type, chr, start, stop, strand, desc]
end

csv.each do |line|
	s = line.split(',')
  id = s[0].gsub('"','')
  if line.start_with?('""')
		tmp_out << %w(ID geneName bioType chromosome start stop strand baseMean log2FoldChange lfcSE pvalue padj fullDescription).join(',') << "\n"
  else
		tmp_out << [id, h[id][0], h[id][1], h[id][2], h[id][3], h[id][4], h[id][5], s[1], s[2], s[3], s[4], s[5].chomp, h[id][6]].join(',') << "\n"
  end    
end
tmp_out.close; csv.close; hash.close
