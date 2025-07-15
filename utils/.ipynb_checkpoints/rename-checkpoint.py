
# coding: utf-8

# In[ ]:


import os
import sys

def changeFileName(dirContainingFile2Rename,extOfFile2Rename,oldString,newString): 
    for f in os.scandir(dirContainingFile2Rename):
        if f.is_dir():
            print('folder found in', f.path ,'!')
        if f.is_file():
            print(os.path.splitext(f.name)[1].lower())
            if os.path.splitext(f.name)[1].lower() in extOfFile2Rename:
                print("oldString is" ,oldString,"newString is" ,newString,'\nf.name is',f.name)
                old_name = f.name
                new_name = old_name.replace(oldString,newString)
                print("rename into == ",new_name,' ==')
                os.rename(f.path, os.path.join(dirContainingFile2Rename,new_name))
            else:
                print('non-'+ext,'files found in', f.path ,'!')
                
                
if (len(sys.argv)-1)!=4:
    print('four arg expected:\n',
          '1st is dirContainingFile2Rename\n',
          '2nd is the extension of files to be renamed\n',
          '3rd is the old string need to be replaced in file name\n',
          '4th is the string above old string will be replaced into\n')
else:
    changeFileName(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

