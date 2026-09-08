import '@jbrowse/react-linear-genome-view2/esm/workerPolyfill'
import { initializeWorker } from '@jbrowse/product-core'
import { enableStaticRendering } from 'mobx-react'
import corePlugins from '@jbrowse/react-linear-genome-view2/esm/corePlugins'
enableStaticRendering(true)

class MyPlugin {
  name = 'MyPlugin'
  install() {
    console.log('myPlugin')
  }
  configure() {}
}
initializeWorker([...corePlugins, MyPlugin], {
  fetchESM: url => import(url),
})

export default function doNothing() {
}
