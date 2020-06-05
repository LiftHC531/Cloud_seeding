# -*- coding: utf-8 -*-
import urllib.request as req
from bs4 import BeautifulSoup

#Python-Cookie
#Cheome: DevTools->application->Cookies(find). 
#Network->Request Headers->cookie 

def getData(url):
    #Build a request, Argument the info. of headers
    request = req.Request(url, headers={
              'cookie':'over18=1',
              'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/81.0.4044.138 Safari/537.36',
    })

    with req.urlopen(request) as response:
         data = response.read().decode("utf-8")

    #print(data)
    root = BeautifulSoup(data, "html.parser")
    #titles = root.find("li", class_="last"); #print(titles.a)
    #titles = root.find_all("li", class_="last") #get a list
    titles = root.find_all("div", class_="title") #get a list
    for title in titles:
        if title.a != None:
           print(title.a.string)
    nextLink =  root.find("a", string="‹ 上頁")
    return nextLink["href"]

#PageURL = "http://radar.atm.ncu.edu.tw/web/#research-four-0-1"
PageURL = "https://www.ptt.cc/bbs/Gossiping/index.html"; #print(str(PageURL)[0:18])
count = 0
while count < 5:
      PageURL = str(PageURL)[0:18]+getData(PageURL)
      print(PageURL)
      count += 1
