import os,sys,optparse,csv

if __name__ == "__main__":
    parser = optparse.OptionParser(usage='usage: %prog runlog.csv [opts] ', version='%prog 1.0')
    parser.add_option('-m', '--runMin'   , type='int'       , default=0 , help='minimum run number')
    parser.add_option('-M', '--runMax'   , type='int'       , default=99999 , help='maximum run number')
    parser.add_option('-o', '--output'   , type='string'    , default="filtered_runlog.csv" , help='output cleaned run log')
    parser.add_option('-H', '--includeHeader'   , action = 'store_true', default=False, help='include the header row (typically when remaking the full CSV)')

    (options, args) = parser.parse_args()

    inputRunlog = args[0]
    outputRunlog = open(options.output,'w')

    with open(inputRunlog,"r") as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
        # This skips the first row (header) of the CSV file.
        if not options.includeHeader:
            next(csvreader)
        for row in list(csvreader):
            if options.includeHeader and row[0]=="run_number":
                outputRunlog.write(' , '.join(row[:-6])+'\n')
            else:
                run=int(row[0])
                if run<options.runMin or run>options.runMax or any('NULL' in field for field in row[:-7]): continue
                outputRunlog.write(' , '.join(row[:-6])+'\n')

    outputRunlog.close()
    print("Filtered runlog in ",options.output)
