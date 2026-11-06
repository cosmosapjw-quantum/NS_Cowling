function Code(el)
 if FORMAT:match('latex') and #el.text > 18 and not el.text:find('[\128-\255]') and not el.text:find('|',1,true) then
 return pandoc.RawInline('latex','\\path|' .. el.text .. '|') end
end