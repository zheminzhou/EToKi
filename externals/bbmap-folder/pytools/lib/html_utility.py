def html_tag(tag, txt, attrs=None):
    html = '<%s'
    if attrs:
        for name in attrs:
            html += ' %s="%s"' % (name, attrs[name])
    html += '>%s</%s>\n'
    html = html % (tag, txt, tag)

    return html

def html_th(header, attrs=None):
    html = ''
    for d in header:
        html += html_tag('th', d, attrs=attrs)
    return html_tag('tr', html)

def html_tr(data):
    html = ''
    for d in data:
        html += html_tag('td', d)
    return html_tag('tr', html)

def html_link(href, txt):
    html = '<a href="%s" target="_blank">%s</a>' % (href, txt)
    return html