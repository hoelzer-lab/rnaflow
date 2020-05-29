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
  id = s[0]; name = s[1]; type = s[2].chomp
  h[id] = [name, type]
end

csv.each do |line|
	s = line.split(',')
  id = s[0].gsub('"','')
  if line.start_with?('""')
		tmp_out << %w(ensemblID geneName bioType baseMean log2FoldChange lfcSE stat pvalue padj).join(',') << "\n"
  else
		tmp_out << [id, h[id][0], h[id][1], s[1], s[2], s[3], s[4], s[5], s[6].chomp].join(',') << "\n"
  end    
end
tmp_out.close; csv.close; hash.close
