# -*- coding: utf-8 -*-
import time
import os
from requests_html import HTMLSession
import urllib
import requests

#Based on pyppeteer (Javescript render)
def get_img(img_urls,model_name):
    for index, img_url in enumerate(img_urls):
        img_element = img_url.split('/');
        print(img_element)
        folder_path = './' + img_element[6][0:10]+'/'
        if (os.path.exists(folder_path) == False):
            os.makedirs(folder_path)
        img_name = model_name+'_'+img_element[7]
        print(img_name)
        urllib.request.urlretrieve(img_url, os.path.join(folder_path, img_name))
        print('第 %d 張' % (index + 1))
        print('output_path:'+folder_path)
        time.sleep(float(index + 1) ** 0.5)

def ncdr_scraper():
    session = HTMLSession()
    url = f"https://watch.ncdr.nat.gov.tw/watch_drought"
    r = session.get(url, headers={
        'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/81.0.4044.138 Safari/537.36',
        'From':'jack661260@gmail.com'
    })
    r.html.render()
    #---NCDR_WRF---
    wrf_urls = r.html.xpath("//div[@style='float:top;width:160px;height:200px;overflow:hidden;']/img/@src")
    #imglink2 = r.html.find("div#float:top;width:160px;height:200px;overflow:hidden; src")
    #---NCDR_MPAS---
    mpas_urls = r.html.xpath("//div[@style='float:top;width:160px;height:230px;overflow:hidden;']/img/@src")
    #months = r.html.search('Python 2 will retire in only {months} months!')['months']
    #print(r.html.html)
    print(wrf_urls)
    #print(imglink2)
    #folder_path.strip()
    get_img(wrf_urls,'WRF')
    get_img(mpas_urls,'MPAS')



ncdr_scraper()
