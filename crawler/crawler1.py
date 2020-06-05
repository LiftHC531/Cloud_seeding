import urllib.request as req
from bs4 import BeautifulSoup

#Chrome:  DevTools->Network->"click"->index.html->Request Headers->user-agent

url = "http://radar.atm.ncu.edu.tw/web/#research-four-0-1"
#Build a request, Argument the info. of headers
request = req.Request(url, headers={
          "User-Agent":"Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/81.0.4044.138 Safari/537.36"})

with req.urlopen(request) as response:
     data = response.read().decode("utf-8")

#print(data)
root = BeautifulSoup(data, "html.parser")
#titles = root.find("li", class_="last"); #print(titles.a)
titles = root.find_all("li", class_="last") #get a list
for title in titles:
    if title.a != None:
       print(title.a.string)
