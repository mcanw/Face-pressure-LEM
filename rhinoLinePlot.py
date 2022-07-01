# Import points from a text file
import csv
import rhinoscriptsyntax as rs
# Imput number of line segments
# number_of_line_segment = rs.GetInteger("Enter number of line segments")
def ImportLines():
    global number_of_line_segment
    number_of_line_segment = 0
    #prompt the user for a file to import
    filter = "CSV file (*.csv)|*.csv|*.txt|All Files (*.*)|*.*||"
    filename = rs.OpenFileName("Open Point File", filter)
    if not filename: return

    with open(filename) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            number_of_line_segment = number_of_line_segment + 1
        
    #read each line from the file
    file = open(filename, "r")
    Startpoint = file.readlines()
    file.close()
 
    #prompt the user for a file to import
    filter = "CSV file (*.csv)|*.csv|All Files (*.*)|*.*||"
    filename = rs.OpenFileName("Open Line File", filter)
    if not filename: return
    
    #read each line from the file
    file = open(filename, "r")
    Endpoint = file.readlines()
    file.close() 
    
    # local helper function    
    def __point_from_string(text):
        items = text.strip("()\n").split(",")
        x = float(items[0])
        y = float(items[1])
        z = float(items[2])
        return x, y, z

    Startpoint = [__point_from_string(line) for line in Startpoint]
    Endpoint = [__point_from_string(line) for line in Endpoint]
    for i in range(0,number_of_line_segment):     
        line = [Startpoint[i],Endpoint[i]]
        lineID = rs.AddLine(line[0],line[1])
    
##########################################################################
# Check to see if this file is being executed as the "main" python
# script instead of being used as a module by some other python script
# This allows us to use the module which ever way we want.
if( __name__ == "__main__" ):
    ImportLines()